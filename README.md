# Pancreatic Cancer Prognosis Using Clinical and Radiomic Data

An educational machine-learning mini-project exploring survival prognosis in pancreatic ductal adenocarcinoma (CPTAC-PDA) using clinical variables and a CT-radiomics feasibility workflow.

> **Research/education only.** This repository is not a diagnostic device, does not provide medical advice, and must not be used to make clinical decisions.

## Current verified status

The project currently contains two separate pieces of work:

1. **Clinical-only survival baseline — executed.** A penalized Cox proportional-hazards model uses age, sex indicator, and tumor-grade group on 143 cleaned clinical cases.
2. **CT radiomics feasibility workflow — executed on three matched patients.** The local extraction output contains three patient rows and 30 numeric radiomics feature columns (plus identifiers/quality-control fields). The RTSTRUCT reference check for C3N-00302 found all 29 referenced CT SOP instances in the recovered CT files, with zero missing references and zero rows in the extraction-error log.

**A combined clinical-plus-radiomics survival model has not been evaluated.** Three matched patients, all with observed events and no censored observations, are far too few for meaningful model evaluation. The repository therefore reports the clinical-only baseline separately and presents radiomics as a feasibility demonstration. Do not describe the combined model as trained or validated.

## Clinical dataset and outcome

The cleaned clinical file is `data/cptac_pda_clinical_survival_clean.csv`:
- 143 usable cases from a source export with 166 rows.
- 115 observed deaths and 28 right-censored observations.
- `survival_days`: days to death for cases recorded as dead; days to last follow-up for cases recorded as alive.
- `event`: 1 for observed death; 0 for right-censoring.
- Cases with unknown/unreported vital status, missing duration, or non-positive duration were excluded.

The model predictors are age at diagnosis in years, a female indicator, and a high-grade indicator (G3/G4 grouped together). See [data/DATA_DICTIONARY.md](data/DATA_DICTIONARY.md).

## Clinical-only baseline results

| Metric | Result |
|---|---:|
| Patients | 143 |
| Observed deaths | 115 |
| Right-censored | 28 |
| Mean 5-fold C-index | 0.567 |
| Fold C-index standard deviation | 0.043 |
| Pooled out-of-fold C-index | 0.574 |

The model uses a penalized Cox proportional-hazards objective and five-fold shuffled cross-validation (seed 42, penalizer 0.1). These are exploratory internal estimates on one dataset with limited predictors—not external validation, proof of clinical usefulness, or a reproduction of the reference paper's clinical-plus-radiomic results.

## Repository structure

- `data/`: cleaned clinical cohort, imaging-overlap/candidate tables, data dictionary, and a radiomics manifest template.
- `src/train_model.py`: clinical-only Cox model and cross-validation.
- `src/extract_radiomics_batch.py`: custom batch feature extraction from locally downloaded CT DICOMs and RTSTRUCTs.
- `src/extract_basic_radiomics.py`: single-case, geometry-matched basic feature extraction.
- `src/verify_radiomics_final.py`: read-only check of saved outputs and CT SOP references for C3N-00302.
- `src/prepare_radiomics_cohort.py`: joins radiomics rows to clinical data and writes a cohort-readiness decision; it does not fit a combined model.
- `src/extract_radiomics.py`: separate PyRadiomics workflow that reads a user-prepared image/mask manifest.
- `src/download_extract_tcia_radiomics.py`: optional TCIA inventory/download/MIRP workflow. Downloading imaging can require several GB and may take a long time.
- `outputs/`: saved model metrics, plots, radiomics CSVs, extraction log, cohort-readiness report, and audit outputs.
- `reports/`: project report PDF, presentation deck, and editable Markdown outlines.
- `.github/workflows/python-check.yml`: syntax-check workflow for selected Python scripts.

Raw CT DICOM downloads, recovered CT folders, NIfTI masks, virtual environments, and other large local imaging files are intentionally excluded from Git. Scripts that process local scans require those files to be obtained separately and placed at the paths expected by the script.

## Run the clinical model (Windows PowerShell)

Use Python 3.10 or newer. From the repository root:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe .\src\train_model.py
```

The script writes `outputs/metrics.json`, `outputs/cross_validation_folds.csv`, `outputs/clinical_cox_coefficients.csv`, and descriptive plots under `outputs/`.

## Radiomics: inspect saved results

The checked-in outputs allow you to inspect the previously generated results without downloading CT data:

```powershell
Get-Content .\outputs\radiomics_cohort_readiness.txt
Import-Csv .\outputs\radiomics_features_cohort.csv | Select-Object patient_id, ct_slices, mask_voxels | Format-Table -AutoSize
```

Current matched rows in the saved feature table:

| Patient ID | CT slices reported by extractor | Mask voxels |
|---|---:|---:|
| C3N-00303 | 161 | 12,833 |
| C3N-00512 | 17 | 9,755 |
| C3N-00302 | 33 | 13,449 |

These are extracted ROI/CT summary features for a small feasibility subset; they are not a validated prognostic radiomics model. The readiness report records 3 matched rows, 3 events, 0 censored observations, and recommends **not** fitting a combined survival model.

### Optional local imaging workflow

The radiomics scripts need local imaging and annotation files. Install the optional dependencies with:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-radiomics.txt
```

Then follow the comments/docstrings in the specific script you intend to run. Do not run the download/extraction workflow during a short presentation: it may download several GB and consume substantial memory. The read-only verifier expects the original `tciaDownload/` folder and the recovered `outputs/selected_ct_series/C3N-00302/` directory to exist locally; those large imaging inputs are not stored in this repository.

## Reproducibility and limitations

- Do not infer that imaging-manifest ID overlap guarantees a usable image/mask pair.
- CT/RTSTRUCT matching and referenced-slice availability are quality-control checks, not proof of complete 3D contour validity.
- The matched radiomics subset is only three patients, all events; it cannot support meaningful combined survival-model evaluation.
- Feature definitions and acquisition differences require additional quality control before any future model comparison.
- Do not report a combined-model score, clinical validation, or clinical benefit from this project.

## Data sources

- [TCIA CPTAC-PDA collection](https://www.cancerimagingarchive.net/collection/cptac-pda/)
- [TCIA CPTAC-PDA tumor annotations](https://www.cancerimagingarchive.net/analysis-result/cptac-pda-tumor-annotations/)
- PDC CPTAC-PDA clinical export used for the cleaned cohort (see data dictionary).
- [Related MultiOmicPDA research repository](https://github.com/FraBer95/MultiOmicPDA)

## License and data handling

No license is specified in this repository yet; check the source dataset's terms before redistributing data. Only share patient-level or imaging-derived data when permitted by the applicable data-use terms.