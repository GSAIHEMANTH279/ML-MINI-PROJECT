"""Download a small matched CPTAC-PDA CT/RTSTRUCT subset and extract MIRP radiomics."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METADATA_URL = "https://www.cancerimagingarchive.net/wp-content/uploads/Metadata_Report_CPTAC_PDA_2025_10_20-1.csv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical", default=str(ROOT / "data" / "cptac_pda_clinical_survival_clean.csv"))
    ap.add_argument("--limit", type=int, default=30, help="Maximum patients; DICOM downloads may require several GB")
    ap.add_argument("--inventory-only", action="store_true", help="Only write candidate list; do not download scans")
    ap.add_argument("--download-only", action="store_true", help="Download DICOM series but skip extraction")
    ap.add_argument("--output", default=str(ROOT / "data" / "radiomics_features.csv"))
    args = ap.parse_args()
    if args.limit < 1:
        raise SystemExit("--limit must be at least 1")
    try:
        from tcia_utils import nbia
    except ImportError:
        raise SystemExit("Install tcia-utils: python -m pip install -r requirements-radiomics.txt")

    print("Reading official TCIA annotation metadata...")
    ann = pd.read_csv(METADATA_URL)
    clinical = pd.read_csv(args.clinical)
    required = {"PatientID", "Annotation Type", "SeriesInstanceUID", "ReferencedSeriesInstanceUID"}
    if not required.issubset(ann.columns):
        raise SystemExit(f"Unexpected TCIA metadata columns; missing {sorted(required-set(ann.columns))}")
    ids = set(clinical["patient_id"].astype(str).str.strip())
    ann["PatientID"] = ann["PatientID"].astype(str).str.strip()
    ann = ann[ann["PatientID"].isin(ids)]
    ann = ann[ann["Annotation Type"].astype(str).str.lower().eq("segmentation")].copy()
    if ann.empty:
        raise SystemExit("No segmentation annotations match the clinical patient IDs.")

    print("Querying TCIA series metadata to keep CT reference images...")
    series = nbia.getSeries("CPTAC-PDA", format="df")
    if not isinstance(series, pd.DataFrame):
        series = pd.DataFrame(series)
    if not {"SeriesInstanceUID", "Modality"}.issubset(series.columns):
        raise SystemExit(f"Unexpected TCIA series columns: {list(series.columns)}")
    modalities = series[["SeriesInstanceUID", "Modality"]].drop_duplicates().copy()
    modalities["SeriesInstanceUID"] = modalities["SeriesInstanceUID"].astype(str)
    ann["ReferencedSeriesInstanceUID"] = ann["ReferencedSeriesInstanceUID"].astype(str)
    ann = ann.merge(modalities.rename(columns={"SeriesInstanceUID":"ReferencedSeriesInstanceUID",
                                               "Modality":"ReferencedModality"}),
                    on="ReferencedSeriesInstanceUID", how="inner")
    ann = ann[ann["ReferencedModality"].astype(str).str.upper().eq("CT")].copy()
    if ann.empty:
        raise SystemExit("No matched CT segmentation rows. Check patient IDs and TCIA metadata.")

    if "ClinicalTrialTimePointID" in ann:
        ann["_tp"] = (~ann["ClinicalTrialTimePointID"].astype(str).str.lower()
                      .str.contains("pre-dose|predose", regex=True)).astype(int)
    else:
        ann["_tp"] = 1
    ann["_date"] = pd.to_datetime(ann["StudyDate"], errors="coerce") if "StudyDate" in ann else pd.NaT
    ann = ann.sort_values(["_tp", "_date", "PatientID"], na_position="last").drop_duplicates("PatientID")
    ann["patient_id"] = ann["PatientID"].astype(str)
    ann["image_series_uid"] = ann["ReferencedSeriesInstanceUID"].astype(str)
    ann["rtstruct_series_uid"] = ann["SeriesInstanceUID"].astype(str)
    ann["image_path"] = ann["image_series_uid"].map(lambda x: str(Path("tciaDownload") / x))
    ann["mask_path"] = ann["rtstruct_series_uid"].map(lambda x: str(Path("tciaDownload") / x / "1-1.dcm"))
    cols = ["patient_id","image_series_uid","rtstruct_series_uid","image_path","mask_path"]
    cols += [c for c in ["StudyDate","ClinicalTrialTimePointID","SeriesDescription",
                          "ReferencedSeriesDescription","ROIVolume","Tracking ID"] if c in ann]
    candidates = ann[cols].reset_index(drop=True)
    inventory = ROOT / "data" / "tcia_matched_annotation_candidates.csv"
    candidates.to_csv(inventory, index=False)
    print(f"Matched clinical patients with annotated CT: {len(candidates)}")
    print(f"Candidate table: {inventory}")
    print(candidates.head(10).to_string(index=False))
    if args.inventory_only:
        print("Inventory-only mode: no scans downloaded.")
        return

    chosen = candidates.head(args.limit)
    uids = sorted(set(chosen["image_series_uid"]) | set(chosen["rtstruct_series_uid"]))
    print(f"Downloading {len(uids)} DICOM series for {len(chosen)} patients. This may take time and several GB.")
    nbia.downloadSeries(uids, input_type="list", format="df")
    if args.download_only:
        print("Download-only mode complete.")
        return

    try:
        from mirp import extract_features
    except ImportError:
        raise SystemExit("Install MIRP: python -m pip install -r requirements-radiomics.txt")
    rows, errors = [], []
    for _, row in chosen.iterrows():
        pid = str(row["patient_id"])
        img, mask = Path(row["image_path"]), Path(row["mask_path"])
        print(f"Extracting {pid}...")
        if not img.is_dir() or not mask.is_file():
            errors.append({"patient_id":pid,"error":"Image series folder or RTSTRUCT file not found"})
            continue
        try:
            result = extract_features(image=str(img), mask=str(mask), by_slice=False,
                                      base_feature_families=["morphology","intensity","texture"])
            frames = [result] if isinstance(result, pd.DataFrame) else [x for x in result if isinstance(x, pd.DataFrame)]
            if not frames:
                raise ValueError("MIRP returned no feature DataFrame")
            numeric = pd.concat(frames, ignore_index=True).select_dtypes(include=[np.number])
            numeric = numeric.replace([np.inf,-np.inf], np.nan)
            values = numeric.median(axis=0, skipna=True)
            row_out = {"patient_id":pid}
            row_out.update({str(k):float(v) for k,v in values.items() if pd.notna(v)})
            if len(row_out) == 1:
                raise ValueError("No numeric radiomics features returned")
            rows.append(row_out)
        except Exception as exc:
            errors.append({"patient_id":pid,"error":f"{type(exc).__name__}: {exc}"})
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f"Saved {len(rows)} patients' radiomics features: {out}")
    if errors:
        err_path = ROOT / "data" / "tcia_radiomics_errors.csv"
        pd.DataFrame(errors).to_csv(err_path, index=False)
        print(f"Extraction errors saved: {err_path}")
    if not rows:
        raise SystemExit("No radiomics features extracted. Inspect the error CSV and verify RTSTRUCT/image pairing.")
    print("Next run: python src/train_combined_model.py")

if __name__ == "__main__":
    main()
