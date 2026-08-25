# Oncology Patient Outcomes Analytics — Breast Cancer (TCGA-BRCA)

## Overview
End-to-end data analytics project analyzing survival outcomes for 1,095 breast cancer patients using real-world clinical data from The Cancer Genome Atlas (TCGA-BRCA). The project covers the full pipeline: data acquisition, database design, statistical analysis, and interactive dashboarding.

## Live Dashboard
[View the interactive dashboard on Tableau Public](https://public.tableau.com/app/profile/fatima.shahid5572/viz/oncology_dashboard/OncologyOverview?publish=yes)

## Key Findings
- Overall mortality rate: **13.79%**
- Mortality rate increases progressively with cancer stage, from ~14% in Stage I to ~75% in Stage IV
- Log-rank test confirms survival differences across stages are statistically significant (p < 0.001)

## Tools & Technologies
- **R** (TCGAbiolinks) — data acquisition from GDC Portal
- **SQL** (SQLite) — relational database design (patients, diagnosis, treatments, outcomes tables)
- **Python** (Pandas, Lifelines, SciPy) — data cleaning, Kaplan-Meier survival analysis, statistical testing
- **Tableau Public** — interactive dashboard and visualization

## Methodology
1. Downloaded clinical data for 1,098 TCGA-BRCA patients via R
2. Designed a normalized SQL schema with four related tables
3. Cleaned and merged data in Python, handling missing values
4. Built Kaplan-Meier survival curves stratified by cancer stage
5. Validated findings with a log-rank statistical test
6. Created an interactive Tableau dashboard summarizing key metrics

## Files
- `TCGA_BRCA_DE_analysis.R` — data acquisition script
- `oncology_data_for_tableau.csv` — cleaned dataset
- `oncology_dashboard.twb` — Tableau workbook
- `km_survival_curve.png` — survival analysis visualization

## Author
Fatima Shahid — MSc Cancer Biology and Pharmacology, Istinye University
