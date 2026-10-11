# Pancreatic Cancer Prognosis Using Clinical and Radiomic Data

## Current verified status
This repository contains a run and evaluated **clinical-only survival-analysis baseline** using the uploaded PDC CPTAC-PDA clinical export. The radiomic part is **not yet implemented in the trained model**, because no radiomic feature table was included in the provided files. The TCIA imaging manifest has 110 unique imaging IDs; 73 of the 143 usable clinical survival records have matching IDs in that manifest. ID overlap alone does not equal radiomic data availability.

## Data and outcome
The raw PDC clinical export had 166 rows. We used `survival_days`: `Days to Death` for `Dead`; `Days to Last Follow Up` for `Alive`. `event` is 1 for `Dead` and 0 for `Alive` (right-censored). Unknown/unreported vital status, missing duration, and non-positive duration were excluded.

The cleaned file contains 143 cases (115 deaths, 28 right-censored). This is an exploratory research dataset, not for clinical decisions.

## Model
Penalized Cox proportional hazards model with three clinical predictors: age at diagnosis (years), sex (female indicator), and tumor grade (G3/G4 grouped as high grade).

Five-fold shuffled cross-validation (seed 42) yielded mean fold C-index **0.567** (SD 0.043) and pooled out-of-fold C-index **0.574**. These are internal exploratory estimates on a small cohort with limited predictors. They are not external validation and do not reproduce the Stanford paper's clinical-plus-radiomic model or its reported results.

## Run on Windows
1. Install Python 3.10+.
2. Open a terminal in this folder.
3. `py -m venv .venv`
4. `.venv\\Scripts\\activate`
5. `python -m pip install -r requirements.txt`
6. `python src/train_model.py`

The script writes metrics, cross-validation fold scores, Cox coefficients and plots to `outputs/`.

## Folder structure
- `data/`: cleaned clinical outcomes, matched-ID list, data dictionary
- `src/train_model.py`: reproducible Cox model and 5-fold evaluation
- `outputs/`: metrics and plots
- `reports/`: report and slide deck

## Required next milestone
Obtain or extract actual radiomic features from tumor-annotated CPTAC-PDA imaging, verify patient-ID matching, and compare clinical-only versus radiomic-only versus combined models with leakage-aware cross-validation. Until that step is completed, describe this as a clinical-only baseline rather than a completed clinical-plus-radiomic model.

## Sources
- TCIA CPTAC-PDA collection: https://www.cancerimagingarchive.net/collection/cptac-pda/
- PDC CPTAC PDA Discovery Study clinical export (downloaded by user)
- Jamalian, Stanford CS229 reference paper (provided by user)
- Related research repository: https://github.com/FraBer95/MultiOmicPDA

## Repository artifacts
- [Cleaned clinical survival cohort](data/cptac_pda_clinical_survival_clean.csv)
- [Clinical/imaging ID overlap table](data/cptac_pda_clinical_imaging_overlap.csv) — ID overlap only, not a radiomics table
- [Model metrics](outputs/metrics.json)
- [Data audit and limitations](outputs/data_audit.json)
- [Cross-validation folds](outputs/cross_validation_folds.csv)
- [Cox coefficients and hazard ratios](outputs/clinical_cox_coefficients.csv)
- [C-index visual summary](outputs/clinical_only_c_index.svg)
- [Report outline](reports/report_outline.md)
- [Presentation outline](reports/presentation_outline.md)
- [Python validation workflow](.github/workflows/python-check.yml)

The PDF and PowerPoint versions are currently available in the accompanying project ZIP, not as files in this repository. The Markdown report and presentation outlines are published here.

## Radiomics pipeline (requires real images and masks)
The repository now includes:
- `src/extract_radiomics.py`: extracts PyRadiomics features from each CT volume and corresponding segmentation mask listed in `data/radiomics_manifest.csv`.
- `data/radiomics_manifest_template.csv`: header template for patient IDs and local image/mask paths.
- `src/train_combined_model.py`: compares clinical-only, radiomics-only and combined Cox models after matching IDs; radiomic feature selection is performed within training folds.

Install the optional extraction dependencies with `python -m pip install pyradiomics SimpleITK`. Copy the manifest template to `data/radiomics_manifest.csv`, fill it with actual local paths to images and validated tumor masks, and run:
```bash
python src/extract_radiomics.py
python src/train_combined_model.py
```
The extraction and combined model have not been run here because the raw image volumes and tumor segmentation masks have not been downloaded and validated. Do not report combined-model metrics until the pipeline has run on real extracted features.

- [PyRadiomics extraction script](src/extract_radiomics.py)
- [Clinical/radiomic comparison script](src/train_combined_model.py)
- [Radiomics manifest template](data/radiomics_manifest_template.csv)
- [Optional radiomics dependencies](requirements-radiomics.txt)

- TCIA CPTAC-PDA tumor annotations (official analysis-results entry): https://www.cancerimagingarchive.net/annotated-cptac-datasets-available-on-tcia/
- Related multi-omic radiomics implementation: https://github.com/FraBer95/MultiOmicPDA
