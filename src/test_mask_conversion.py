
from pathlib import Path
from collections import Counter, defaultdict
import shutil
import pydicom
import SimpleITK as sitk
from rt_utils import RTStructBuilder

ROOT = Path("tciaDownload")
OUT = Path("outputs")
OUT.mkdir(exist_ok=True)

# Find RTSTRUCT
rt_files = []
for p in ROOT.rglob("*.dcm"):
    try:
        d = pydicom.dcmread(str(p), stop_before_pixels=True, force=True)
        if getattr(d, "Modality", "") == "RTSTRUCT":
            rt_files.append(p)
    except Exception:
        pass

if not rt_files:
    raise RuntimeError("No RTSTRUCT found.")

rt_path = rt_files[0]
rt_ds = pydicom.dcmread(str(rt_path), force=True)

series_uid = str(
    rt_ds.ReferencedFrameOfReferenceSequence[0]
    .RTReferencedStudySequence[0]
    .RTReferencedSeriesSequence[0].SeriesInstanceUID
)

# SOP Instance UIDs referenced by the RTSTRUCT contours
referenced_sops = set()
for roi in getattr(rt_ds, "ROIContourSequence", []):
    for contour in getattr(roi, "ContourSequence", []):
        for image_ref in getattr(contour, "ContourImageSequence", []):
            referenced_sops.add(
                str(image_ref.ReferencedSOPInstanceUID)
            )

print("RTSTRUCT:", rt_path)
print("Referenced CT series:", series_uid)
print("Referenced contour images:", len(referenced_sops))

# Index CT files in the referenced series by acquisition
groups = defaultdict(list)
for p in ROOT.rglob("*.dcm"):
    try:
        d = pydicom.dcmread(str(p), stop_before_pixels=True, force=True)
        if (
            getattr(d, "Modality", "") == "CT"
            and str(getattr(d, "SeriesInstanceUID", "")) == series_uid
        ):
            acquisition = str(getattr(d, "AcquisitionNumber", "unknown"))
            groups[acquisition].append((p, d))
    except Exception:
        pass

if not groups:
    raise RuntimeError("Referenced CT series not found.")

# Select the acquisition containing the most contour-referenced images
scores = {
    acq: sum(
        str(getattr(d, "SOPInstanceUID", "")) in referenced_sops
        for _, d in items
    )
    for acq, items in groups.items()
}
print("Acquisition scores:", scores)

best_acq = max(scores, key=scores.get)
if scores[best_acq] == 0:
    raise RuntimeError(
        "No CT acquisition matches the RTSTRUCT contour image references. "
        "Do not generate a mask until this is resolved."
    )

selected = groups[best_acq]
print("Selected acquisition:", best_acq)
print("Selected CT files:", len(selected))

# Keep only the selected acquisition in a separate directory
filtered = OUT / "selected_ct_series"
if filtered.exists():
    shutil.rmtree(filtered)
filtered.mkdir(parents=True)

for p, _ in selected:
    shutil.copy2(p, filtered / p.name)

# Read the selected CT volume
reader = sitk.ImageSeriesReader()
names = reader.GetGDCMSeriesFileNames(str(filtered))
if not names:
    raise RuntimeError("SimpleITK could not read the selected CT series.")
reader.SetFileNames(names)
ct = reader.Execute()

print("CT size:", ct.GetSize())
print("CT spacing:", ct.GetSpacing())

# Create the RTSTRUCT mask against the selected CT acquisition
rtstruct = RTStructBuilder.create_from(
    dicom_series_path=str(filtered),
    rt_struct_path=str(rt_path),
)

roi_names = rtstruct.get_roi_names()
print("Available ROIs:", roi_names)

roi_name = next(
    (name for name in roi_names if "PANCREAS" in name.upper()),
    None,
)
if roi_name is None:
    raise RuntimeError("No pancreas ROI found.")

mask = rtstruct.get_roi_mask_by_name(roi_name)
mask_array = mask.transpose(2, 0, 1).astype("uint8")

ct_array = sitk.GetArrayFromImage(ct)
print("CT array shape (z,y,x):", ct_array.shape)
print("Mask array shape (z,y,x):", mask_array.shape)

if ct_array.shape != mask_array.shape:
    raise RuntimeError(
        "CT and mask dimensions differ. Stop: alignment is not verified."
    )

mask_image = sitk.GetImageFromArray(mask_array)
mask_image.CopyInformation(ct)

mask_path = OUT / "pancreas_mask_aligned.nii.gz"
sitk.WriteImage(mask_image, str(mask_path))

print("ROI:", roi_name)
print("Mask voxels:", int(mask_array.sum()))
print("Saved:", mask_path)
print("Mask spacing:", mask_image.GetSpacing())
print("Mask origin:", mask_image.GetOrigin())
print("Mask direction:", mask_image.GetDirection())
print("Dimension and geometry metadata copied from CT.")
