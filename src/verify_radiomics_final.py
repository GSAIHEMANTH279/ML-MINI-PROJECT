"""Final verification for the pancreatic cancer prognosis mini-project.

Run from the project root:
    python src/verify_radiomics_final.py

This is read-only: it does not change cohort outputs.
"""
from pathlib import Path
import csv
import sys

import pydicom

ROOT = Path(".")
DICOM_ROOT = ROOT / "tciaDownload"
RECOVERED = ROOT / "outputs" / "selected_ct_series" / "C3N-00302"
FEATURES = ROOT / "outputs" / "radiomics_features_cohort.csv"
ERRORS = ROOT / "outputs" / "radiomics_extraction_errors.csv"
MERGED = ROOT / "outputs" / "clinical_radiomics_survival_cohort.csv"
READINESS = ROOT / "outputs" / "radiomics_cohort_readiness.txt"

TARGET = "C3N-00302"


def read_dicom(path):
    return pydicom.dcmread(str(path), stop_before_pixels=True, force=True)


def get_references(rt):
    series_uids, sop_uids = set(), set()
    try:
        for frame in rt.ReferencedFrameOfReferenceSequence:
            for study in frame.RTReferencedStudySequence:
                for series in study.RTReferencedSeriesSequence:
                    series_uids.add(str(series.SeriesInstanceUID))
                    for item in getattr(series, "ContourImageSequence", []):
                        uid = getattr(item, "ReferencedSOPInstanceUID", None)
                        if uid:
                            sop_uids.add(str(uid))
    except Exception:
        pass
    try:
        for roi in rt.ROIContourSequence:
            for contour in getattr(roi, "ContourSequence", []):
                for item in getattr(contour, "ContourImageSequence", []):
                    uid = getattr(item, "ReferencedSOPInstanceUID", None)
                    if uid:
                        sop_uids.add(str(uid))
    except Exception:
        pass
    return series_uids, sop_uids


def load_csv(path):
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    print("=== FINAL PROJECT VERIFICATION (read-only) ===")
    required = [FEATURES, ERRORS, MERGED, READINESS]
    for path in required:
        print(f"{'EXISTS' if path.exists() else 'MISSING'}: {path}")
    if not FEATURES.exists():
        print("Cannot verify feature rows because the feature CSV is missing.")
        return 2

    rows = load_csv(FEATURES)
    ids = [str(r.get("patient_id", "")).strip() for r in rows]
    print(f"\nFeature rows: {len(rows)}")
    print(f"Patient IDs: {', '.join(ids) if ids else '(none)'}")
    for row in rows:
        print(
            f"  {row.get('patient_id')}: ct_slices={row.get('ct_slices')}, "
            f"mask_voxels={row.get('mask_voxels')}, "
            f"features/columns={len(row)}"
        )

    errors = load_csv(ERRORS)
    print(f"\nExtraction error rows: {len(errors)}")
    for err in errors:
        print(f"  {err.get('patient_id')}: {err.get('error')}")

    # Verify that all RTSTRUCT-referenced SOP instances for C3N-00302 are
    # present in the recovered series folder.
    if not RECOVERED.exists():
        print(f"\nMASK-COVERAGE CHECK: recovered folder missing: {RECOVERED}")
        return 2

    ct_by_sop = {}
    for path in RECOVERED.rglob("*.dcm"):
        try:
            ds = read_dicom(path)
            if str(getattr(ds, "Modality", "")).upper() != "CT":
                continue
            sop = str(getattr(ds, "SOPInstanceUID", ""))
            if sop:
                ct_by_sop[sop] = (
                    path,
                    str(getattr(ds, "SeriesInstanceUID", "")),
                    float(ds.ImagePositionPatient[2]) if hasattr(ds, "ImagePositionPatient") else None,
                )
        except Exception:
            continue
    print(f"\nRecovered CT images readable: {len(ct_by_sop)}")

    rt_candidates = []
    if DICOM_ROOT.exists():
        for path in DICOM_ROOT.rglob("*.dcm"):
            try:
                ds = read_dicom(path)
                if str(getattr(ds, "Modality", "")).upper() == "RTSTRUCT" and str(getattr(ds, "PatientID", "")) == TARGET:
                    rt_candidates.append((path, ds))
            except Exception:
                continue

    if not rt_candidates:
        print(f"MASK-COVERAGE CHECK: could not locate RTSTRUCT with PatientID {TARGET} under {DICOM_ROOT}.")
        print("Feature extraction succeeded, but complete ROI coverage cannot be confirmed by this check.")
        return 1

    rt_path, rt_ds = rt_candidates[0]
    series_uids, ref_sops = get_references(rt_ds)
    found_refs = ref_sops & set(ct_by_sop)
    missing_refs = ref_sops - set(ct_by_sop)
    print(f"RTSTRUCT: {rt_path}")
    print(f"Referenced CT Series UIDs: {', '.join(sorted(series_uids)) or '(none)'}")
    print(f"Unique referenced CT SOP UIDs: {len(ref_sops)}")
    print(f"Referenced SOP UIDs present in recovered CT: {len(found_refs)}")
    print(f"Referenced SOP UIDs missing from recovered CT: {len(missing_refs)}")

    if ref_sops:
        positions = [ct_by_sop[sop][2] for sop in found_refs if ct_by_sop[sop][2] is not None]
        if positions:
            print(f"Referenced contour z-position span (mm): {min(positions):.3f} to {max(positions):.3f}")
        if missing_refs:
            print("First missing referenced SOP UIDs:")
            for uid in sorted(missing_refs)[:10]:
                print(f"  {uid}")
    else:
        print("WARNING: RTSTRUCT exposes no slice-level CT SOP references; full mask coverage cannot be established from UID matching.")

    print("\n=== INTERPRETATION ===")
    if TARGET not in ids:
        print(f"FAIL: {TARGET} is not in the feature CSV.")
        return 1
    if missing_refs:
        print("WARNING: some referenced CT instances are absent; do not claim complete RTSTRUCT slice coverage.")
        return 1
    if not ref_sops:
        print("INCOMPLETE: no slice-level SOP references available to verify complete mask coverage.")
        return 1
    if errors:
        print("WARNING: extraction error log is not empty.")
        return 1
    if not MERGED.exists() or not READINESS.exists():
        print("WARNING: one or more cohort output files are missing.")
        return 1

    print("PASS: all referenced CT SOP instances were found and required output files exist.")
    print("NOTE: this verifies referenced-slice availability, not clinical validity or full 3D contour completeness.")
    print("NOTE: three matched patients are insufficient to evaluate a combined survival model.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
