# Fetal Cardiac Radiomics Pipeline

Analytical code and de-identified dataset for the study **"An exploration of ultrasound-based radiomics combined with machine learning in the diagnosis of fetal ventricular disproportion caused by cardiovascular malformations"**.

This pipeline implements feature selection and repeated stratified cross-validation for:
1. Three-class one-vs-all classification (Control / GroupA / GroupB)
2. Binary classification of GroupA vs GroupB using radiomics, clinical, and combined feature sets
3. Sensitivity analyses restricted to (a) left heart obstructive lesions and (b) cases with postnatal echocardiographic confirmation

---

## Repository structure

```
.
├── analysis_pipeline.R       # Main analysis script
├── data/
│   ├── radiomics_dataset.xlsx  # De-identified dataset
│   └── ICC_select.csv          # Features passing ICC reproducibility threshold
├── output/                   # Generated automatically when the script is run
└── README.md
```

---

## Requirements

- R version 4.3.3 or later
- The following R packages (installed automatically if missing):
  `readxl`, `dplyr`, `pROC`, `glmnet`, `ranger`, `xgboost`, `tidyr`

---


## Data dictionary

| Variable | Description |
|---|---|
| `No` | Case identifier (re-assigned for sharing; no link to original medical records) |
| `Group` | Diagnostic group: "0"(Contro), "1" (GroupA：ventricular disproportion without cardiovascular malformations), "2" (GroupB：ventricular disproportion with cardiovascular malformations) |
| `Disease` | Specific cardiovascular malformation for GroupB cases |
| `PostnatalEcho` | 1 = case with institutional postnatal echocardiographic confirmation; 0 = case followed up at referring institutions |
| `Instrument` | Ultrasound machine used "0" (E10) or "1" (E8))|
| `GA` | Gestational age at examination (weeks) |
| `LV_GLS` | Left ventricular global longitudinal strain (%) |
| `MV.annulus` | Mitral valve annulus diameter (mm) |
| `RV.LV` | Right-to-left ventricular mid-transverse diameter ratio |
| Columns from `AD` onward | Radiomic features extracted from the right ventricular free wall, interventricular septum, and left ventricular lateral wall |

---

## Output files

Key outputs (all CSV) include:
- `3class_OVA_summary.csv`, `binary_AvsB_summary.csv` — model performance metrics
- `*_feature_frequency.csv` — within-fold feature selection frequencies
- `*_DeLong.csv` — pairwise model comparisons
- `*_calibration_summary.csv`, `*_calibration_bins.csv` — calibration metrics and binned data for plotting
- `*_OOF_predictions.csv` — out-of-fold predicted probabilities
- `*_repeat_preds.csv` — per-repetition OOF predictions for ROC confidence band plotting
- `comparison_*_main_vs_SA.csv`, `echo_SA_vs_main_comparison.csv` — primary vs sensitivity analysis comparisons


## Reproducibility note

Due to the inherent stochasticity of cross-validation, exact numerical reproduction of the reported metrics may show minor deviations (typically <0.02 in AUC). The overall conclusions, model rankings, and statistical comparisons remain consistent.