
from pathlib import Path
from collections import Counter, defaultdict
import pydicom

uid = (
    "1.3.6.1.4.1.14519.5.2.1."
    "3320.3273.228476454931090952243913193877"
)

by_z = defaultdict(list)

for path in Path("tciaDownload").rglob("*.dcm"):
    try:
        ds = pydicom.dcmread(
            str(path), stop_before_pixels=True, force=True
        )
        if (
            getattr(ds, "Modality", "") == "CT"
            and str(getattr(ds, "SeriesInstanceUID", "")) == uid
            and hasattr(ds, "ImagePositionPatient")
        ):
            z = round(float(ds.ImagePositionPatient[2]), 4)
            by_z[z].append({
                "file": str(path),
                "instance": getattr(ds, "InstanceNumber", "?"),
                "sop": str(getattr(ds, "SOPInstanceUID", "?")),
            })
    except Exception:
        pass

positions = sorted(by_z)
gaps = [
    round(positions[i + 1] - positions[i], 4)
    for i in range(len(positions) - 1)
]

print("Files:", sum(len(v) for v in by_z.values()))
print("Unique positions:", len(positions))
print("Duplicate-position count:", sum(len(v) > 1 for v in by_z.values()))
print("Common gaps (mm):", Counter(gaps).most_common(10))

print("\nDuplicate positions:")
for z, items in by_z.items():
    if len(items) > 1:
        print(z, [(x["instance"], x["file"]) for x in items])

print("\nGaps larger than 5 mm:")
for i, gap in enumerate(gaps):
    if gap > 5:
        print(positions[i], "to", positions[i + 1], "=", gap, "mm")
