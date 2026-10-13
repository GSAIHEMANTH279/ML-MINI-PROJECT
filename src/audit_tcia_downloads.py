"""Audit locally downloaded CPTAC-PDA DICOMs without downloading anything."""
from pathlib import Path
from collections import defaultdict
import csv
import pydicom

ROOT = Path('tciaDownload')
OUT = Path('outputs') / 'dicom_download_audit.csv'
if not ROOT.exists():
    raise SystemExit(f'Cannot find {ROOT.resolve()}. Run from the project root folder.')

ct_series = defaultdict(lambda: {'patients': set(), 'acquisitions': set(), 'files': 0, 'paths': set()})
rtstructs = []
errors = 0
for p in ROOT.rglob('*.dcm'):
    try:
        ds = pydicom.dcmread(str(p), stop_before_pixels=True, force=True)
        modality = str(getattr(ds, 'Modality', '')).upper()
        patient = str(getattr(ds, 'PatientID', 'UNKNOWN'))
        uid = str(getattr(ds, 'SeriesInstanceUID', ''))
        if modality == 'CT' and uid:
            key = (patient, uid, str(getattr(ds, 'AcquisitionNumber', 'NA')))
            row = ct_series[key]
            row['patients'].add(patient); row['acquisitions'].add(str(getattr(ds, 'AcquisitionNumber', 'NA')))
            row['files'] += 1; row['paths'].add(str(p.parent))
        elif modality == 'RTSTRUCT':
            ref_uids = set()
            try:
                for frame in ds.ReferencedFrameOfReferenceSequence:
                    for study in frame.RTReferencedStudySequence:
                        for series in study.RTReferencedSeriesSequence:
                            ref_uids.add(str(series.SeriesInstanceUID))
            except Exception:
                pass
            rois = []
            try:
                rois = [str(x.ROIName) for x in ds.StructureSetROISequence]
            except Exception:
                pass
            rtstructs.append({'path': str(p), 'patient_id': patient, 'series_uid': uid,
                              'referenced_ct_uids': ';'.join(sorted(ref_uids)),
                              'roi_names': ';'.join(rois)})
    except Exception:
        errors += 1

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open('w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['record_type','patient_id','series_uid','acquisition_number','file_count','folder','referenced_ct_uids','roi_names'])
    for (patient, uid, acq), r in sorted(ct_series.items()):
        w.writerow(['CT', patient, uid, acq, r['files'], ';'.join(sorted(r['paths'])), '', ''])
    for r in rtstructs:
        w.writerow(['RTSTRUCT', r['patient_id'], r['series_uid'], '', '', r['path'], r['referenced_ct_uids'], r['roi_names']])

print('LOCAL DICOM AUDIT COMPLETE')
print('Root:', ROOT.resolve())
print('CT patient/series/acquisition groups:', len(ct_series))
print('Distinct CT PatientIDs:', len({k[0] for k in ct_series}))
print('RTSTRUCT files:', len(rtstructs))
print('Distinct RTSTRUCT PatientIDs:', len({r['patient_id'] for r in rtstructs}))
print('DICOM read errors skipped:', errors)
for r in rtstructs:
    print('\nRTSTRUCT patient:', r['patient_id'])
    print('  file:', r['path'])
    print('  referenced CT series:', r['referenced_ct_uids'] or '(not found in reference sequence)')
    print('  ROIs:', r['roi_names'] or '(none listed)')
print('\nSaved audit CSV:', OUT)
if len(rtstructs) < 2:
    print('NOTE: Fewer than 2 local RTSTRUCTs found. A multi-patient radiomics model cannot be trained from the current local download alone.')
