# Dataset and outcome definitions

Source: user-exported PDC clinical tables for CPTAC PDA Discovery Study and TCIA CPTAC-PDA imaging metadata manifest.

## Files
- `cptac_pda_clinical_survival_clean.csv`: 143 clinical cases with valid positive survival duration and known vital status.
- `cptac_pda_clinical_imaging_overlap.csv`: IDs for the 73 cleaned clinical cases also appearing in the TCIA imaging metadata manifest.
- `../outputs/data_audit.json`: audit counts, model configuration, and limitations.

## Target
- `survival_days`: `Days to Death` for cases with source vital status `Dead`; `Days to Last Follow Up` for cases with source vital status `Alive`.
- `event`: `1` for `Dead` (observed event); `0` for `Alive` (right-censored).
- Rows with missing/unreported vital status, missing survival duration, or non-positive duration are excluded.

## Predictors in the executed model
- `age_years`: Age at Diagnosis converted from days to years (divide by 365.25).
- `female`: derived in code from source `Sex`.
- `high_grade`: derived in code as G3/G4 versus other grades.

## Radiomics status
No radiomic feature values were included in the clinical export. Although 73 clinical cases match IDs in the imaging metadata manifest, this is not enough to claim radiomic integration. The executed results are clinical-only. Do not label them as the full clinical-plus-radiomic model.
