# Dataset and outcome definitions

Source: user-exported PDC clinical tables for the CPTAC-PDA Discovery Study and TCIA CPTAC-PDA imaging/annotation metadata.

## Files
- `cptac_pda_clinical_survival_clean.csv`: 143 clinical cases with valid positive survival duration and known vital status.
- `cptac_pda_clinical_imaging_overlap.csv`: cleaned clinical IDs that also appear in the TCIA imaging metadata manifest.
- `tcia_matched_annotation_candidates.csv`: candidate imaging/annotation matches identified by the optional TCIA workflow.
- `../outputs/data_audit.json`: clinical data audit, model configuration, and limitations.
- `../outputs/radiomics_features_cohort.csv`: extracted CT/ROI summary-feature table for three matched patients.
- `../outputs/radiomics_cohort_readiness.txt`: join counts and decision on whether the matched cohort is suitable for combined modeling.

## Target
- `survival_days`: `Days to Death` for source vital status `Dead`; `Days to Last Follow Up` for `Alive`.
- `event`: 1 for `Dead` (observed event); 0 for `Alive` (right-censored).
- Rows with missing/unreported vital status, missing duration, or non-positive duration are excluded.

## Predictors in the executed clinical model
- `age_years`: Age at Diagnosis converted from days to years (divide by 365.25).
- `female`: derived in code from source `Sex`.
- `high_grade`: derived in code as G3/G4 versus other grades.

## Radiomics status
The saved radiomics feature table contains three patient rows and 30 detected radiomics feature columns, excluding identifiers and quality-control columns. All three matched patients have event=1 and none are right-censored. The readiness script therefore advises against fitting/reporting a combined survival model from this cohort. The clinical-only model results must be reported separately.

The verification script's CT reference check confirms availability of the RTSTRUCT-referenced CT SOP instances it checks; it does not prove full 3D contour completeness or clinical validity.

## Interpretation and data use
These files support educational and exploratory research only. Do not claim clinical validation. Share patient-level clinical or imaging-derived data only when the source license and data-use terms permit redistribution.