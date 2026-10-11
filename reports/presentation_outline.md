# Presentation Outline

## Slide 1 — Title
Pancreatic Cancer Prognosis Using Clinical and Radiomic Data.

## Slide 2 — Objective
Investigate clinical predictors of survival in pancreatic ductal adenocarcinoma and prepare for comparison with imaging-derived radiomic features.

## Slide 3 — Data
166 clinical records in the source export; 143 usable survival records after cleaning; 115 observed deaths and 28 right-censored cases.

## Slide 4 — Imaging Linkage
110 unique patient IDs in the TCIA manifest; 73 IDs overlap the cleaned clinical cohort. Matching IDs do not guarantee usable radiomic feature values.

## Slide 5 — Method
Penalized Cox proportional hazards model using age, sex and tumor grade; five-fold shuffled cross-validation with seed 42.

## Slide 6 — Results
Mean fold C-index 0.567 (SD 0.043); pooled out-of-fold C-index 0.574. These are exploratory internal estimates.

## Slide 7 — Limitations
The trained model is clinical-only. Radiomic features and tumor annotations are not yet integrated. Results are not external validation.

## Slide 8 — Next Steps
Obtain real radiomic features for matched patients and compare clinical-only, radiomic-only and combined models.