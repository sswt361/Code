# Project Name: TLS-Related Gene Analysis

## Project Introduction
This project aims to analyze the expression levels of TLS (Tertiary Lymphoid Structures)-related genes and assess their differences between healthy and tumor tissues through various analytical methods. Additionally, it explores the correlation with clinical features.

## Code Description

### Code 7
- **Filename**: `code7.kt` 
- **Purpose**: Extract expression levels of TLS-related genes.
- **Input**: Gene expression dataset.
- **Output**: Expression levels of TLS-related genes.

### Code 12
- **Filename**: `code12.kt`
- **Purpose**: Perform differential analysis between healthy tissue and tumor tissue.
- **Input**: Gene expression data from healthy and tumor tissues.
- **Output**: Statistical results of differentially expressed genes.

### Code 14
- **Filename**: `code14.kt`
- **Purpose**: Conduct ssGSEA (Single Sample Gene Set Enrichment Analysis) on TCGA data.
- **Input**: TCGA gene expression data.
- **Output**: ssGSEA analysis results.

### Code 15
- **Filename**: `code15.kt` 
- **Purpose**: Perform TLS score survival analysis.
- **Input**: TLS score data and patient survival data.
- **Output**: Survival analysis results (e.g., Kaplan-Meier curves).

### Code 16
- **Filename**: `code16.kt`
- **Purpose**: Conduct ssGSEA analysis on GEO data.
- **Input**: GEO gene expression data.
- **Output**: ssGSEA analysis results.

### Code 18
- **Filename**: `code18.kt` 
- **Purpose**: Perform clinical differential analysis based on ssGSEA analysis results.
- **Input**: ssGSEA analysis results and clinical feature data.
- **Output**: Results of clinical differential analysis.

## Project Structure
TLS_Analysis/
├── code7.py
├── code12.py
├── code14.py
├── code15.py
├── code16.py
├── code18.py
├── data/
│   ├── TCGA_data.csv
│   ├── GEO_data.csv
│   └── clinical_data.csv
├── results/
│   ├── expression_levels.csv
│   ├── differential_analysis.csv
│   ├── ssGSEA_results_TCGA.csv
│   ├── ssGSEA_results_GEO.csv
│   └── survival_analysis_results.csv
└── README.md
## Usage
1. Place the data files in the `data/` folder.
2. Run the corresponding code files to generate results.
3. Results will be stored in the `results/` folder.

## Dependencies and Environment
- R
