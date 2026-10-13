# Mini-Project Presentation Outline

## Slide 1 — Title
Pancreatic Cancer Prognosis Using Clinical and Radiomic Data.

## Slide 2 — Objective
Explore clinical survival predictors and demonstrate a CT-radiomics extraction workflow for pancreatic cancer imaging.

## Slide 3 — Clinical data
The source export contained 166 rows; 143 cases had valid survival duration and known vital status after cleaning. There were 115 observed deaths and 28 right-censored observations.

## Slide 4 — Clinical model
Penalized Cox proportional-hazards model with age, female indicator, and high-grade indicator. Five-fold shuffled cross-validation, seed 42, penalizer 0.1.

## Slide 5 — Clinical-only results
Mean fold C-index: 0.567 (SD 0.043). Pooled out-of-fold C-index: 0.574. These are exploratory internal estimates, not external validation.

## Slide 6 — CT radiomics feasibility
Saved feature table: three matched patients and 30 detected radiomics feature columns. The checked C3N-00302 RTSTRUCT had 29 referenced CT SOP instances present in the recovered CT folder, with zero missing references.

## Slide 7 — Limitations
The three matched radiomics cases all have observed events and none are censored. This is too small to evaluate a combined survival model. Referenced-slice availability does not prove complete 3D segmentation validity. The executed survival model is clinical-only.

## Slide 8 — Conclusion and next steps
The clinical-only baseline runs reproducibly and the radiomics extraction is demonstrated as a feasibility workflow. Future work needs a larger, quality-controlled matched cohort before comparing clinical-only, radiomics-only, and combined models.

## Presentation note
Do not claim the combined model was trained, validated, or achieved a performance score.