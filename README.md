# Pancreatic Cancer Prognosis Using Clinical and Radiomic Data

## Current verified status
This repository now contains a **run and evaluated clinical-only survival-analysis baseline** using the uploaded PDC CPTAC-PDA clinical export. The radiomic part is **not yet implemented in the trained model**, because no radiomic feature table was included in the provided files. The TCIA imaging manifest has 110 unique imaging IDs; 73 of the 143 usable clinical survival records have matching IDs in that manifest. ID overlap alone does not equal radiomic data availability.

## Data and outcome
The raw PDC clinical export had 166 rows. We used:
- `survival_days`: `Days to Death` for `Dead`; `Days to Last Follow Up` for `Alive`.
- `event`: 1 for `Dead`; 0 for `Alive` (right-censored).
- Excluded unknown/unreported vital status, missing duration, and non-positive duration.

The cleaned file contains 143 cases (115 deaths, 28 right-censored). This is an exploratory research dataset, not for clinical decisions.

## Model
Penalized Cox proportional hazards model with three clinical predictors:
- Age at diagnosis (years)
- Sex (female indicator)
- Tumor grade (G3/G4 grouped as high grade)

Five-fold shuffled cross-validation (seed 42) yielded:
- Mean fold C-index: **0.567** (SD 0.043)
- Pooled out-of-fold C-index: **0.574**

These are internal exploratory estimates on a small cohort with limited predictors. They are not external validation and do not reproduce the Stanford paper's clinical-plus-radiomic model or its reported results.

## Run on Windows
1. Install Python 3.10+.
2. Open terminal in this folder.
3. `py -m venv .venv`
4. `.venv\\Scripts\\activate`
5. `python -m pip install -r requirements.txt`
6. `python src/train_model.py`

The script writes metrics, cross-validation fold scores, Cox coefficients and plots to `outputs/`.

## Folder structure
- `data/`: cleaned clinical outcomes, matched-ID list, data dictionary
- `src/train_model.py`: reproducible Cox model and 5-fold evaluation
- `outputs/`: actual metrics and plots from the executed run
- `reports/`: report and slide deck

## Required next milestone for the full project title
Obtain or extract actual radiomic features from tumor-annotated CPTAC-PDA imaging, verify patient-ID matching, and then compare clinical-only versus radiomic-only versus combined models with leakage-aware cross-validation. Until that step is completed, describe this as a clinical-only baseline rather than a completed clinical-plus-radiomic model.

## Sources
- TCIA CPTAC-PDA collection: https://www.cancerimagingarchive.net/collection/cptac-pda/
- PDC CPTAC PDA Discovery Study clinical export (downloaded by user)
- Jamalian, Stanford CS229 reference paper (provided by user)
- MultiOmicPDA related research repository: https://github.com/FraBer95/MultiOmicPDA

## Radiomics pipeline (requires real images and masks)

- `src/extract_radiomics.py` extracts PyRadiomics features from each image volume and segmentation mask listed in `data/radiomics_manifest.csv`.
- `data/radiomics_manifest_template.csv` is a header template for patient IDs and local image/mask paths.
- `src/train_combined_model.py` compares clinical-only, radiomics-only, and combined Cox models after matching IDs. Radiomic feature selection is performed within training folds.
- `requirements-radiomics.txt` lists optional extraction dependencies.

Install extraction dependencies with `python -m pip install -r requirements-radiomics.txt`. Copy the manifest template to `data/radiomics_manifest.csv`, fill it with real local image and validated tumor-mask paths, and run:

```bash
python src/extract_radiomics.py
python src/train_combined_model.py
```

The extraction and combined model have not been run because the raw image volumes and tumor segmentation masks have not been downloaded and validated. Do not report combined-model metrics until the pipeline has run on real extracted features.

## Current repository status
The executed result is a clinical-only Cox baseline: 143 cases, 115 observed deaths, 28 right-censored cases, mean five-fold C-index 0.567 (SD 0.043), pooled out-of-fold C-index 0.574. The report and slide deck are in `reports/`. The GitHub repository includes Markdown report/presentation content; the PDF and PowerPoint are included in the project ZIP.
