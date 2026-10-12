
from pathlib import Path
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt

ct_folder = Path("outputs/selected_ct_series")
mask_path = Path("outputs/pancreas_mask_aligned.nii.gz")

reader = sitk.ImageSeriesReader()
files = reader.GetGDCMSeriesFileNames(str(ct_folder))
reader.SetFileNames(files)
ct = reader.Execute()

mask = sitk.ReadImage(str(mask_path))

ct_arr = sitk.GetArrayFromImage(ct)
mask_arr = sitk.GetArrayFromImage(mask)

if ct_arr.shape != mask_arr.shape:
    raise ValueError(
        f"Shape mismatch: CT {ct_arr.shape}, mask {mask_arr.shape}"
    )

# Display the middle slice containing pancreas mask voxels
slices = np.where(mask_arr.sum(axis=(1, 2)) > 0)[0]
if len(slices) == 0:
    raise ValueError("Mask is empty.")

z = int(slices[len(slices) // 2])

plt.figure(figsize=(8, 8))
plt.imshow(ct_arr[z], cmap="gray", vmin=-150, vmax=250)
overlay = np.ma.masked_where(mask_arr[z] == 0, mask_arr[z])
plt.imshow(overlay, cmap="autumn", alpha=0.55)
plt.title(f"CT with pancreas mask — slice {z}")
plt.axis("off")
plt.tight_layout()
plt.savefig("outputs/pancreas_mask_overlay.png", dpi=160)
plt.show()

print("Saved overlay: outputs/pancreas_mask_overlay.png")
print("Slices containing mask:", len(slices))
