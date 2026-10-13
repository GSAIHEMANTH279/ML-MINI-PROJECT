# Pancreatic Cancer Prognosis — Demo and Submission

## What was actually completed
- Clinical-only Cox proportional hazards baseline on 143 patients (115 observed deaths; 28 censored).
- Five-fold CV mean concordance index 0.567; pooled out-of-fold concordance index 0.574.
- Radiomics extraction feasibility demonstrated for three matched cases, with 30 radiomics feature columns and no logged extraction errors.
- **No validated combined clinical-plus-radiomics survival model**: three matched cases, all events, are insufficient for reliable evaluation.

## Live demo (Windows PowerShell)
Open the project root, activate your existing virtual environment, then run:

```powershell
.\.venv-radiomics\Scripts\Activate.ps1
python src/train_model.py
```

Inspect saved outputs without re-downloading large CT data:

```powershell
Get-Content outputs\RESULTS_SUMMARY.txt
Import-Csv outputs\radiomics_features_cohort.csv | Select-Object -First 3 | Format-Table
Get-Content outputs\radiomics_cohort_readiness.txt
```

If you use a different virtual environment, activate that environment instead. Do not rerun bulk TCIA downloads during the live demo. The radiomics extraction scripts depend on large local CT/RTSTRUCT files not committed to GitHub.

## Explain the workflow
1. Clinical data cleaning: survival duration and event indicators; 143 eligible records.
2. Cox model: age, female indicator and high-grade indicator; five-fold cross-validation.
3. Concordance index: ability to rank survival risk, not accuracy or probability of cure.
4. Radiomics: CT/RTSTRUCT alignment, feature extraction, 3-case feasibility cohort.
5. Limitation: combined model could not be reliably trained or validated with only three matched cases.

## Final submission verification
- Upload the final **two-page PDF report** to `reports/Mini_Project_Report.pdf`.
- Upload the final **PowerPoint presentation** to `reports/Mini_Project_Presentation.pptx`.
- Confirm that both files open after uploading; binary files cannot be represented by Markdown outlines alone.
- If the assignment requires a private repository, change repository visibility to private and explicitly invite faculty/TAs before submission. The repository was public at last check.
- Verify faculty access, submission link, and live demonstration locally.

Official sources: [TCIA CPTAC-PDA](https://www.cancerimagingarchive.net/collection/cptac-pda/) and [TCIA tumor annotations](https://www.cancerimagingarchive.net/analysis-result/cptac-pda-tumor-annotations/).
