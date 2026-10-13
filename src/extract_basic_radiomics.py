"""Extract simple research radiomics from the exact CT acquisition matching the mask.
Run from project root: python src/extract_basic_radiomics.py
Requires: pydicom, SimpleITK, numpy, scipy.
"""
from pathlib import Path
from collections import defaultdict
import csv
import numpy as np
import SimpleITK as sitk
from scipy import ndimage
import pydicom

ROOT = Path('.')
DICOM_ROOT = ROOT / 'tciaDownload'
MASK_PATH = ROOT / 'outputs' / 'pancreas_mask_aligned.nii.gz'
OUT_PATH = ROOT / 'outputs' / 'radiomics_features_single_case.csv'

def close_tuple(a, b, tol=1e-3):
    return len(a) == len(b) and all(abs(float(x)-float(y)) <= tol for x, y in zip(a, b))

def read_ct_acquisition_groups(root):
    """Group individual CT files by SeriesInstanceUID AND AcquisitionNumber.
    This avoids mixing separate acquisitions that share a Series UID in a download.
    """
    groups = defaultdict(list)
    for path in root.rglob('*.dcm'):
        try:
            ds = pydicom.dcmread(str(path), stop_before_pixels=True, force=True)
            if getattr(ds, 'Modality', '') != 'CT' or not hasattr(ds, 'ImagePositionPatient'):
                continue
            uid = str(getattr(ds, 'SeriesInstanceUID', ''))
            acq = str(getattr(ds, 'AcquisitionNumber', 'NA'))
            iop = getattr(ds, 'ImageOrientationPatient', None)
            ipp = getattr(ds, 'ImagePositionPatient', None)
            if not uid or ipp is None or iop is None:
                continue
            row = np.asarray(iop[:3], dtype=float)
            col = np.asarray(iop[3:6], dtype=float)
            normal = np.cross(row, col)
            pos = float(np.dot(np.asarray(ipp, dtype=float), normal))
            groups[(uid, acq)].append((pos, str(path), ds))
        except Exception:
            continue
    candidates = []
    for (uid, acq), items in groups.items():
        # Keep one file for each physical slice position; report groups separately.
        items.sort(key=lambda x: x[0])
        unique = []
        seen = set()
        for item in items:
            pos_key = round(item[0], 4)
            if pos_key not in seen:
                seen.add(pos_key)
                unique.append(item)
        if len(unique) < 10:
            continue
        names = [x[1] for x in unique]
        try:
            reader = sitk.ImageSeriesReader()
            reader.SetFileNames(names)
            image = reader.Execute()
            candidates.append({'uid': uid, 'acq': acq, 'files': names,
                               'image': image, 'folder': str(Path(names[0]).parent)})
        except Exception as exc:
            print(f'Skipping UID={uid}, acquisition={acq}: {exc}')
    return candidates

def main():
    if not MASK_PATH.exists():
        raise FileNotFoundError(f'Mask not found: {MASK_PATH}')
    if not DICOM_ROOT.exists():
        raise FileNotFoundError(f'DICOM folder not found: {DICOM_ROOT}')
    mask_img = sitk.ReadImage(str(MASK_PATH))
    mask = sitk.GetArrayFromImage(mask_img) > 0
    if not mask.any():
        raise ValueError('Mask has zero voxels.')
    print('Mask shape (z,y,x):', mask.shape, '| voxels:', int(mask.sum()))
    print('Grouping CT DICOM files by series UID + acquisition number...')
    candidates = read_ct_acquisition_groups(DICOM_ROOT)
    print('CT acquisition groups readable:', len(candidates))
    matches = []
    for c in candidates:
        im = c['image']
        shape = tuple(reversed(im.GetSize()))
        same = (shape == tuple(mask.shape)
                and close_tuple(im.GetSpacing(), mask_img.GetSpacing(), 1e-3)
                and close_tuple(im.GetOrigin(), mask_img.GetOrigin(), 1e-2)
                and close_tuple(im.GetDirection(), mask_img.GetDirection(), 1e-3))
        print(f"UID={c['uid']} acquisition={c['acq']} slices={shape[0]} size={im.GetSize()} origin={im.GetOrigin()} spacing={im.GetSpacing()} {'MATCH' if same else ''}")
        if same:
            matches.append(c)
    if len(matches) != 1:
        raise RuntimeError(f'Expected exactly one matching CT acquisition; found {len(matches)}. The candidate summary above shows the mismatch. Do not force-extract features from a nonmatching volume.')
    chosen = matches[0]
    ct = sitk.GetArrayFromImage(chosen['image']).astype(np.float32)
    spacing_xyz = chosen['image'].GetSpacing()
    spacing_zyx = (spacing_xyz[2], spacing_xyz[1], spacing_xyz[0])
    vox = ct[mask]
    vox = vox[np.isfinite(vox)]
    if not vox.size:
        raise ValueError('No finite CT intensities inside mask.')
    features = {
        'ct_series_uid': chosen['uid'], 'acquisition_number': chosen['acq'],
        'ct_folder': chosen['folder'], 'mask_path': str(MASK_PATH),
        'mask_voxels': int(mask.sum()), 'voxel_volume_mm3': float(np.prod(spacing_xyz)),
        'mask_volume_mm3': float(mask.sum() * np.prod(spacing_xyz)),
        'intensity_mean': float(np.mean(vox)), 'intensity_median': float(np.median(vox)),
        'intensity_std': float(np.std(vox, ddof=1)) if vox.size > 1 else 0.0,
        'intensity_min': float(np.min(vox)), 'intensity_max': float(np.max(vox)),
        'intensity_range': float(np.ptp(vox)), 'intensity_q10': float(np.percentile(vox, 10)),
        'intensity_q25': float(np.percentile(vox, 25)), 'intensity_q75': float(np.percentile(vox, 75)),
        'intensity_q90': float(np.percentile(vox, 90)),
        'intensity_iqr': float(np.percentile(vox, 75)-np.percentile(vox, 25)),
        'intensity_skewness': float(np.mean(((vox-vox.mean())/(vox.std()+1e-12))**3)),
        'intensity_kurtosis_excess': float(np.mean(((vox-vox.mean())/(vox.std()+1e-12))**4)-3),
    }
    labels, ncomp = ndimage.label(mask)
    counts = np.bincount(labels.ravel())[1:]
    features['connected_components_3d'] = int(ncomp)
    features['largest_component_voxels'] = int(counts.max()) if counts.size else 0
    coords = np.argwhere(mask)
    ext = coords.max(axis=0)-coords.min(axis=0)+1
    features['bbox_extent_z_mm'] = float(ext[0]*spacing_zyx[0])
    features['bbox_extent_y_mm'] = float(ext[1]*spacing_zyx[1])
    features['bbox_extent_x_mm'] = float(ext[2]*spacing_zyx[2])
    features['mask_slices'] = int(np.count_nonzero(mask.reshape(mask.shape[0], -1).sum(axis=1)))
    for name, axis in [('z',0),('y',1),('x',2)]:
        s1=[slice(None)]*3; s2=[slice(None)]*3
        s1[axis]=slice(0,-1); s2[axis]=slice(1,None)
        s1=tuple(s1); s2=tuple(s2)
        valid=mask[s1] & mask[s2] & np.isfinite(ct[s1]) & np.isfinite(ct[s2])
        diffs=ct[s2][valid]-ct[s1][valid]
        features[f'{name}_neighbor_pair_count']=int(diffs.size)
        features[f'{name}_neighbor_abs_diff_mean']=float(np.mean(np.abs(diffs))) if diffs.size else float('nan')
        features[f'{name}_neighbor_diff_std']=float(np.std(diffs)) if diffs.size else float('nan')
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(features.keys())); w.writeheader(); w.writerow(features)
    print('\nFEATURE EXTRACTION SUCCEEDED')
    print('CT series UID:', chosen['uid'], '| acquisition:', chosen['acq'])
    print('Feature count:', len(features)); print('Saved:', OUT_PATH)
    print('Research/demo features only; one case, not a cohort model.')

if __name__ == '__main__':
    main()
