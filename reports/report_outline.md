# Mini-Project Report Outline

## Problem
Pancreatic cancer prognosis using clinical and radiomic data. The completed modeling result is a clinical-only survival baseline; CT radiomics was extracted as a feasibility demonstration but is not integrated into a validated combined model.

## Clinical dataset
The source clinical export contained 166 rows. After excluding unknown/unreported vital status, missing survival duration, and non-positive duration, 143 usable survival cases remained: 115 observed deaths and 28 right-censored cases.

## Method
A penalized Cox proportional-hazards model uses age at diagnosis in years, a female indicator, and a high-grade indicator (G3/G4 grouped). Evaluation used five-fold shuffled cross-validation with random seed 42 and penalizer 0.1. Preprocessing was fit within each training fold.

## Clinical-only results
- Mean fold C-index: 0.567
- Fold C-index standard deviation: 0.043
- Pooled out-of-fold C-index: 0.574

These are exploratory internal estimates, not external validation or evidence of clinical utility.

## Radiomics feasibility
The saved cohort radiomics table contains three matched patient rows and 30 detected radiomics feature columns. The readiness report records three events and zero censored observations among matched rows and explicitly recommends against fitting/reporting a combined survival model from this subset. The verification script found all 29 CT SOP instances referenced by the checked C3N-00302 RTSTRUCT in the recovered CT files; this checks referenced-slice availability, not complete 3D contour validity.

## Limitations
The combined clinical-plus-radiomics model has not been evaluated. Three matched cases are insufficient for meaningful survival-model evaluation. CT acquisition and mask geometry require careful matching, and extracted features are a feasibility output only.

## Sources
- https://www.cancerimagingarchive.net/collection/cptac-pda/
- https://www.cancerimagingarchive.net/analysis-result/cptac-pda-tumor-annotations/
- https://github.com/FraBer95/MultiOmicPDA

## Reproducibility
The repository includes a Python syntax-check workflow at `.github/workflows/python-check.yml`. Run the clinical baseline using the README instructions. Radiomics scripts require local imaging files that are intentionally not committed.

## Data governance
Redistribute clinical and imaging-derived patient-level data only when the source license and applicable data-use terms permit public sharing.