# Mini-Project Report

## Problem
Pancreatic cancer prognosis using clinical and radiomic data. Current experiment is a clinical-only baseline; radiomics is not integrated yet.

## Dataset
166 source clinical records; 143 usable survival cases after excluding missing status/duration and non-positive duration. 115 observed deaths and 28 right-censored cases. 73 patient IDs overlap the TCIA imaging manifest, but that does not establish availability of radiomic features.

## Method
Penalized Cox proportional hazards model; predictors age in years, female indicator, and high-grade indicator. Five-fold shuffled cross-validation, random seed 42, with preprocessing fitted within each training fold.

## Results
Mean fold C-index 0.567 (SD 0.043); pooled out-of-fold C-index 0.574. These are exploratory internal estimates, not clinical validation.

## Limitations
Radiomic features and tumor annotations have not been integrated. Do not claim reproduction of the reference paper's combined model or metrics. Next step is obtain/derive real radiomic features for matched IDs and compare clinical-only, radiomic-only, and combined models.

## Sources
- https://www.cancerimagingarchive.net/collection/cptac-pda/
- https://github.com/FraBer95/MultiOmicPDA
