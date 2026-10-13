"""Join extracted CT radiomics to clinical survival data and report whether modeling is feasible.

Run from the project root:
    python src/prepare_radiomics_cohort.py

Expected inputs:
    outputs/radiomics_features_cohort.csv
    data/cptac_pda_clinical_survival_clean.csv

This script does not invent missing data or fit a survival model when the sample is too small.
"""
from pathlib import Path
import csv
import sys

ROOT = Path(".")
RAD_PATH = ROOT / "outputs" / "radiomics_features_cohort.csv"
CLINICAL_CANDIDATES = [
    ROOT / "data" / "cptac_pda_clinical_survival_clean.csv",
    ROOT / "data" / "clinical_survival_clean.csv",
    ROOT / "data" / "clinical_clean.csv",
]
OUT_DIR = ROOT / "outputs"
MERGED_PATH = OUT_DIR / "clinical_radiomics_survival_cohort.csv"
REPORT_PATH = OUT_DIR / "radiomics_cohort_readiness.txt"


def read_csv(path):
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not RAD_PATH.exists():
        sys.exit(f"Missing {RAD_PATH}. Run the batch radiomics extractor first.")

    clinical_path = next((p for p in CLINICAL_CANDIDATES if p.exists()), None)
    if clinical_path is None:
        sys.exit(
            "Could not find cleaned clinical survival CSV. Expected one of:\n"
            + "\n".join(str(p) for p in CLINICAL_CANDIDATES)
        )

    radiomics = read_csv(RAD_PATH)
    clinical = read_csv(clinical_path)
    if not radiomics or not clinical:
        sys.exit("Radiomics or clinical CSV has no data rows.")

    if "patient_id" not in radiomics[0] or "patient_id" not in clinical[0]:
        sys.exit("Both CSVs must contain a patient_id column.")

    clinical_by_id = {}
    for row in clinical:
        pid = (row.get("patient_id") or "").strip()
        if pid:
            clinical_by_id[pid] = row

    feature_cols = [
        key for key in radiomics[0]
        if key not in {"patient_id", "ct_series_uid", "acquisition_number",
                       "rtstruct_path", "roi_name", "ct_slices", "mask_voxels"}
    ]
    # Keep identifiers and quality-control fields as well as clinical outcome columns.
    merged = []
    unmatched = []
    for rrow in radiomics:
        pid = (rrow.get("patient_id") or "").strip()
        crow = clinical_by_id.get(pid)
        if crow is None:
            unmatched.append(pid)
            continue
        row = dict(crow)
        for key, value in rrow.items():
            if key == "patient_id":
                continue
            row[f"radiomics_{key}"] = value
        merged.append(row)

    if merged:
        fieldnames = list(merged[0].keys())
        with MERGED_PATH.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(merged)

    event_count = sum(str(r.get("event", "")).strip() == "1" for r in merged)
    censored_count = sum(str(r.get("event", "")).strip() == "0" for r in merged)
    report = [
        "RADIOMICS + CLINICAL COHORT READINESS",
        f"Clinical file: {clinical_path}",
        f"Radiomics rows: {len(radiomics)}",
        f"Clinical rows: {len(clinical)}",
        f"Matched rows: {len(merged)}",
        f"Radiomics patient IDs absent from clinical data: {', '.join(unmatched) if unmatched else 'none'}",
        f"Events among matched rows: {event_count}",
        f"Censored observations among matched rows: {censored_count}",
        f"Radiomics feature columns detected: {len(feature_cols)}",
        "",
    ]
    if len(merged) < 20:
        report += [
            "DECISION: Do not fit or report a combined survival model from this cohort.",
            "The matched sample is too small for a meaningful model evaluation.",
            "Use the radiomics extraction as a feasibility demonstration and report the clinical-only baseline separately.",
        ]
    else:
        report += [
            "DECISION: The cohort has at least 20 matched rows, but this alone does not establish model adequacy.",
            "Before modeling, check event count, missingness, feature leakage, and use nested/train-fold-only feature selection.",
        ]

    REPORT_PATH.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    if merged:
        print(f"\nSaved merged cohort: {MERGED_PATH}")
    print(f"Saved readiness report: {REPORT_PATH}")


if __name__ == "__main__":
    main()
