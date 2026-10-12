"""Batch-extract simple 3D pancreatic-mask radiomics from locally downloaded CPTAC-PDA DICOMs.
Run from the project root: python src/extract_radiomics_batch.py
Requires: pydicom, rt-utils, SimpleITK, numpy, scipy.
Only extracts when the RTSTRUCT can be matched to one unambiguous CT acquisition.
"""
from pathlib import Path
from collections import defaultdict
import csv, shutil, tempfile
import numpy as np
import pydicom
import SimpleITK as sitk
from scipy import ndimage
from rt_utils import RTStructBuilder

ROOT = Path('.')
DICOM_ROOT = ROOT / 'tciaDownload'
RECOVERED_C3N302 = ROOT / 'outputs' / 'selected_ct_series' / 'C3N-00302'
OUT_DIR = ROOT / 'outputs'
OUT_CSV = OUT_DIR / 'radiomics_features_cohort.csv'
ERROR_CSV = OUT_DIR / 'radiomics_extraction_errors.csv'


def dicom_read(path, pixels=False):
    return pydicom.dcmread(str(path), stop_before_pixels=not pixels, force=True)


def referenced_ct_info(ds):
    """Return referenced CT SeriesInstanceUID and referenced CT SOPInstanceUIDs."""
    series_uids, sop_uids = set(), set()
    try:
        for frame in ds.ReferencedFrameOfReferenceSequence:
            for study in frame.RTReferencedStudySequence:
                for series in study.RTReferencedSeriesSequence:
                    series_uids.add(str(series.SeriesInstanceUID))
                    for inst in getattr(series, 'ContourImageSequence', []):
                        if getattr(inst, 'ReferencedSOPInstanceUID', None):
                            sop_uids.add(str(inst.ReferencedSOPInstanceUID))
    except Exception:
        pass
    # SOP references are usually stored under ROI contour image sequences.
    try:
        for roi in ds.ROIContourSequence:
            for contour in getattr(roi, 'ContourSequence', []):
                for image in getattr(contour, 'ContourImageSequence', []):
                    if getattr(image, 'ReferencedSOPInstanceUID', None):
                        sop_uids.add(str(image.ReferencedSOPInstanceUID))
    except Exception:
        pass
    return series_uids, sop_uids


def feature_row(patient_id, rt_path, ct_uid, acq, ct_img, mask_img, mask):
    ct = sitk.GetArrayFromImage(ct_img).astype(np.float32)
    vox = ct[mask]
    vox = vox[np.isfinite(vox)]
    if not vox.size:
        raise ValueError('No finite CT voxels inside pancreas mask')
    sp = ct_img.GetSpacing()  # x,y,z
    row = {
        'patient_id': patient_id, 'ct_series_uid': ct_uid, 'acquisition_number': acq,
        'rtstruct_path': str(rt_path), 'mask_voxels': int(mask.sum()),
        'voxel_volume_mm3': float(np.prod(sp)), 'mask_volume_mm3': float(mask.sum()*np.prod(sp)),
        'intensity_mean': float(np.mean(vox)), 'intensity_median': float(np.median(vox)),
        'intensity_std': float(np.std(vox, ddof=1)) if vox.size > 1 else 0.0,
        'intensity_min': float(np.min(vox)), 'intensity_max': float(np.max(vox)),
        'intensity_range': float(np.ptp(vox)), 'intensity_q10': float(np.percentile(vox,10)),
        'intensity_q25': float(np.percentile(vox,25)), 'intensity_q75': float(np.percentile(vox,75)),
        'intensity_q90': float(np.percentile(vox,90)),
        'intensity_iqr': float(np.percentile(vox,75)-np.percentile(vox,25)),
        'intensity_skewness': float(np.mean(((vox-vox.mean())/(vox.std()+1e-12))**3)),
        'intensity_kurtosis_excess': float(np.mean(((vox-vox.mean())/(vox.std()+1e-12))**4)-3),
    }
    _, ncomp = ndimage.label(mask)
    labels, counts = np.unique(ndimage.label(mask)[0], return_counts=True)
    counts = counts[labels != 0]
    row['connected_components_3d'] = int(ncomp)
    row['largest_component_voxels'] = int(counts.max()) if counts.size else 0
    coords = np.argwhere(mask); ext = coords.max(axis=0)-coords.min(axis=0)+1
    row['bbox_extent_z_mm'] = float(ext[0]*sp[2]); row['bbox_extent_y_mm'] = float(ext[1]*sp[1]); row['bbox_extent_x_mm'] = float(ext[2]*sp[0])
    row['mask_slices'] = int(np.count_nonzero(mask.reshape(mask.shape[0], -1).sum(axis=1)))
    for name, axis in [('z',0),('y',1),('x',2)]:
        s1=[slice(None)]*3; s2=[slice(None)]*3; s1[axis]=slice(0,-1); s2[axis]=slice(1,None)
        s1=tuple(s1); s2=tuple(s2); valid=mask[s1] & mask[s2] & np.isfinite(ct[s1]) & np.isfinite(ct[s2])
        diffs=ct[s2][valid]-ct[s1][valid]
        row[f'{name}_neighbor_pair_count']=int(diffs.size)
        row[f'{name}_neighbor_abs_diff_mean']=float(np.mean(np.abs(diffs))) if diffs.size else float('nan')
        row[f'{name}_neighbor_diff_std']=float(np.std(diffs)) if diffs.size else float('nan')
    return row


def main():
    if not DICOM_ROOT.exists():
        raise SystemExit(f'Missing {DICOM_ROOT}. Run from project root.')
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ct_groups = defaultdict(list); rt_paths=[]
    # Scan original downloads plus the recovered CT series for C3N-00302.
    # RTSTRUCT files are still discovered from DICOM_ROOT; the recovery folder
    # contains CT DICOMs only.
    scan_paths = list(DICOM_ROOT.rglob('*.dcm'))
    if RECOVERED_C3N302.exists():
        scan_paths.extend(RECOVERED_C3N302.rglob('*.dcm'))
    else:
        print(f'NOTE: recovered CT folder not found: {RECOVERED_C3N302}')
    seen_paths = set()
    for p in scan_paths:
        resolved = str(p.resolve())
        if resolved in seen_paths:
            continue
        seen_paths.add(resolved)
        try:
            ds=dicom_read(p)
            modality=str(getattr(ds,'Modality','')).upper()
            if modality=='CT' and hasattr(ds,'ImagePositionPatient') and hasattr(ds,'SeriesInstanceUID'):
                key=(str(getattr(ds,'PatientID','UNKNOWN')),str(ds.SeriesInstanceUID),str(getattr(ds,'AcquisitionNumber','NA')))
                ct_groups[key].append(p)
            elif modality=='RTSTRUCT': rt_paths.append(p)
        except Exception:
            continue
    print(f'Found {len(ct_groups)} CT acquisition groups and {len(rt_paths)} RTSTRUCT files.')
    successful=[]; errors=[]; seen_patients=set()
    for rt_path in rt_paths:
        try:
            rtds=dicom_read(rt_path)
            patient=str(getattr(rtds,'PatientID','UNKNOWN'))
            series_uids, sop_uids=referenced_ct_info(rtds)
            if not series_uids:
                raise ValueError('RTSTRUCT has no referenced CT SeriesInstanceUID')
            candidates=[]
            for (pid,uid,acq), files in ct_groups.items():
                if pid != patient or uid not in series_uids: continue
                file_sops=set()
                for f in files:
                    try: file_sops.add(str(dicom_read(f).SOPInstanceUID))
                    except Exception: pass
                overlap=len(file_sops & sop_uids)
                candidates.append((overlap,(pid,uid,acq),files,file_sops))
            if not candidates:
                # Some DICOMs have inconsistent PatientID; series UID remains the authoritative match.
                for (pid,uid,acq),files in ct_groups.items():
                    if uid not in series_uids: continue
                    file_sops=set()
                    for f in files:
                        try: file_sops.add(str(dicom_read(f).SOPInstanceUID))
                        except Exception: pass
                    candidates.append((len(file_sops & sop_uids),(pid,uid,acq),files,file_sops))
            if not candidates: raise ValueError(f'No downloaded CT group matches referenced UID(s): {sorted(series_uids)}')
            if len(candidates) > 1 and not sop_uids:
                raise ValueError('Multiple CT acquisitions share the referenced Series UID, but RTSTRUCT has no slice-level SOP references; refusing ambiguous match')
            candidates.sort(key=lambda x:(x[0],len(x[2])), reverse=True)
            best=candidates[0]
            if len(candidates)>1 and best[0]==candidates[1][0]:
                raise ValueError('Ambiguous CT acquisition: referenced SOP instances do not uniquely identify acquisition')
            if sop_uids and best[0] == 0:
                raise ValueError('No referenced CT SOP instances overlap candidate files; refusing unsafe match')
            _, (ct_patient,ct_uid,acq), files, file_sops=best

            # rt-utils may allocate a float64 mask for every CT slice. For
            # very long series, retain the contiguous span containing the
            # RTSTRUCT-referenced SOP instances, plus two-slice margins.
            if len(files) > 300 and sop_uids:
                positioned = []
                for f in files:
                    try:
                        d = dicom_read(f)
                        z = float(d.ImagePositionPatient[2])
                        sop = str(getattr(d, "SOPInstanceUID", ""))
                        positioned.append((z, f, sop))
                    except Exception:
                        continue
                positioned.sort(key=lambda item: item[0])
                ref_indices = [
                    i for i, (_, _, sop) in enumerate(positioned) if sop in sop_uids
                ]
                if ref_indices:
                    lo = max(0, min(ref_indices) - 2)
                    hi = min(len(positioned), max(ref_indices) + 3)
                    cropped = [item[1] for item in positioned[lo:hi]]
                    if len(cropped) >= 3 and len(cropped) < len(files):
                        print(
                            f"MEMORY-SAFE CT subset for {patient}: {len(cropped)} "
                            f"contiguous slices of {len(files)} total"
                        )
                        files = cropped
                    else:
                        print(
                            f"NOTE {patient}: referenced slice span did not reduce "
                            f"the series; retaining {len(files)} slices."
                        )

            if (ct_patient,ct_uid,acq) in seen_patients:
                print(f'SKIP {patient}: already extracted a case for this patient')
                continue
            # Keep only selected acquisition files in a temporary folder so rt-utils cannot mix acquisitions.
            with tempfile.TemporaryDirectory(prefix='pancreas_ct_') as tmp:
                tmpdir=Path(tmp)
                for i, f in enumerate(files): shutil.copy2(f, tmpdir / f'slice_{i:05d}.dcm')
                # Read CT volume from precisely this acquisition.
                reader=sitk.ImageSeriesReader()
                series_ids=sitk.ImageSeriesReader.GetGDCMSeriesIDs(str(tmpdir))
                if not series_ids: raise ValueError('SimpleITK cannot read selected CT acquisition')
                # Use files from the series ID with the most files; temp contains only this acquisition.
                names=sitk.ImageSeriesReader.GetGDCMSeriesFileNames(str(tmpdir),series_ids[0])
                reader.SetFileNames(list(names)); ct_img=reader.Execute()
                rt=RTStructBuilder.create_from(dicom_series_path=str(tmpdir),rt_struct_path=str(rt_path))
                roi_names=rt.get_roi_names()
                roi=next((n for n in roi_names if 'PANCREAS' in n.upper()),None)
                if roi is None: raise ValueError(f'No pancreas ROI found; available ROIs={roi_names}')
                mask_xyz=rt.get_roi_mask_by_name(roi)
                mask_zyx=np.transpose(mask_xyz,(2,0,1)).astype(np.uint8)
                if tuple(mask_zyx.shape) != tuple(reversed(ct_img.GetSize())):
                    raise ValueError(f'Mask/CT shape mismatch: mask {mask_zyx.shape}, CT {tuple(reversed(ct_img.GetSize()))}')
                mask_img=sitk.GetImageFromArray(mask_zyx); mask_img.CopyInformation(ct_img)
                if not mask_zyx.any(): raise ValueError('Extracted pancreas mask is empty')
                row=feature_row(patient,rt_path,ct_uid,acq,ct_img,mask_img,mask_zyx.astype(bool))
                row['roi_name']=roi; row['ct_slices']=int(ct_img.GetSize()[2])
                successful.append(row); seen_patients.add((ct_patient,ct_uid,acq))
                print(f"SUCCESS {patient}: ROI={roi}, slices={row['ct_slices']}, voxels={row['mask_voxels']}, features={len(row)}")
        except Exception as exc:
            errors.append({'patient_id':str(getattr(locals().get('rtds',None),'PatientID','UNKNOWN')),'rtstruct_path':str(rt_path),'error':str(exc)})
            print(f'FAILED {rt_path.name}: {exc}')
    # Merge successful new rows with existing output rows by patient ID. This
    # prevents a later failed extraction from deleting previously successful cases.
    existing_rows = []
    if OUT_CSV.exists():
        try:
            with OUT_CSV.open('r', newline='', encoding='utf-8-sig') as f:
                existing_rows = list(csv.DictReader(f))
        except Exception as exc:
            print(f'WARNING: could not read previous cohort CSV: {exc}')

    merged_by_patient = {
        str(row.get('patient_id', '')): row
        for row in existing_rows
        if str(row.get('patient_id', '')).strip()
    }
    for row in successful:
        merged_by_patient[str(row['patient_id'])] = row
    merged_rows = list(merged_by_patient.values())

    if merged_rows:
        fields = list(dict.fromkeys(k for row in merged_rows for k in row.keys()))
        with OUT_CSV.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(merged_rows)
        print(f'\nSaved {len(merged_rows)} patient rows (new successes this run: {len(successful)}): {OUT_CSV}')
    else:
        print('\nNo radiomics rows extracted; no output cohort CSV created.')
    with ERROR_CSV.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['patient_id','rtstruct_path','error']); w.writeheader(); w.writerows(errors)
    print(f'Saved extraction log: {ERROR_CSV}')
    print('NOTE: these are simple research/demo features, not a clinically validated radiomics pipeline.')

if __name__=='__main__': main()
