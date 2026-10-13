
from pathlib import Path
from collections import Counter
import pydicom

# CT series referenced by your pancreas annotation
uid = "1.3.6.1.4.1.14519.5.2.1.3320.3273.228476454931090952243913193877"

rows = []
for path in Path("tciaDownload").rglob("*.dcm"):
    try:
        d = pydicom.dcmread(
            str(path), stop_before_pixels=True, force=True
        )
        if (
            getattr(d, "Modality", "") == "CT"
            and str(getattr(d, "SeriesInstanceUID", "")) == uid
            and hasattr(d, "ImagePositionPatient")
        ):
            rows.append((float(d.ImagePositionPatient[2]), str(path)))
    except Exception:
        pass

rows.sort()
positions = [round(z, 4) for z, _ in rows]
gaps = [
    round(positions[i + 1] - positions[i], 4)
    for i in range(len(positions) - 1)
]

print("CT files with positions:", len(rows))
print("Unique slice positions:", len(set(positions)))
print("Duplicate files:", len(rows) - len(set(positions)))
print("Most common slice gaps (mm):", Counter(gaps).most_common(8))
print("Gaps greater than 5 mm:", [
    (positions[i], positions[i + 1], gaps[i])
    for i in range(len(gaps)) if gaps[i] > 5
][:15])
