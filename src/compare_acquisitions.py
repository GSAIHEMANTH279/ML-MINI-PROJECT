from pathlib import Path
from collections import defaultdict
import pydicom

uid = "1.3.6.1.4.1.14519.5.2.1.3320.3273.228476454931090952243913193877"
groups = defaultdict(list)

for p in Path("tciaDownload").rglob("*.dcm"):
    try:
        d = pydicom.dcmread(str(p), stop_before_pixels=True, force=True)
        if (getattr(d, "Modality", "") == "CT"
                and str(getattr(d, "SeriesInstanceUID", "")) == uid
                and hasattr(d, "ImagePositionPatient")):
            key = str(getattr(d, "AcquisitionNumber", "missing"))
            groups[key].append(float(d.ImagePositionPatient[2]))
    except Exception:
        pass

for key, values in sorted(groups.items()):
    unique = sorted(set(round(z, 3) for z in values))
    gaps = [round(unique[i+1] - unique[i], 3)
            for i in range(len(unique)-1)]
    print(f"\nAcquisition {key}")
    print("Files:", len(values))
    print("Unique positions:", len(unique))
    print("Z range:", (unique[0], unique[-1]) if unique else "none")
    print("Common gaps:", __import__("collections").Counter(gaps).most_common(5))
