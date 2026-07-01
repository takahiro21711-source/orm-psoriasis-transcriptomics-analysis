# ORM1/ORM2 Psoriasis Transcriptomics Analysis

## Overview
This project provides a comprehensive bioinformatics pipeline for analyzing **Orosomucoid 1 and 2 (ORM1/ORM2)** gene expression in psoriasis using harmonized transcriptomics data from Zenodo.

### Key Features
- 📊 **Data Integration**: Harmonized psoriasis transcriptomics dataset (39 datasets, 1,677 samples)
- 🔬 **Gene Expression Analysis**: ORM1/ORM2 expression profiling
- 📈 **Differential Expression**: Psoriasis vs healthy controls comparison
- 🧬 **Pathway Analysis**: GO enrichment and KEGG pathway analysis
- 🎯 **Epigenetic Analysis**: TFBS prediction and chromatin accessibility inference
- 📉 **Visualization**: Interactive plots, heatmaps, network diagrams, HTML reports

## Dataset
- **Source**: Zenodo DOI: 10.5281/zenodo.4009497
- **Publication**: Federico et al., Sci Data 7, 153 (2020)
- **Data**: Manually curated and harmonized transcriptomics datasets of psoriasis
- **Samples**: 1,677 samples from 39 different studies
- **Content**: Gene expression matrices with clinical metadata

## Installation

### Requirements
- Python 3.9 or higher
- pip or conda

### Setup
```bash
# Clone the repository
git clone https://github.com/takahiro21711-source/orm-psoriasis-transcriptomics-analysis.git
cd orm-psoriasis-transcriptomics-analysis

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Project Structure
```
orm-psoriasis-transcriptomics-analysis/
├── README.md
├── requirements.txt
├── config.yaml
├── data/
│   └── .gitkeep
├── notebooks/
│   ├── 01_data_download_preprocessing.ipynb
│   ├── 02_orm_expression_analysis.ipynb
│   ├── 03_differential_expression.ipynb
│   ├── 04_pathway_enrichment.ipynb
│   ├── 05_epigenetic_analysis.ipynb
│   └── 06_integrated_visualization.ipynb
├── scripts/
│   ├── download_data.py
│   ├── data_preprocessing.py
│   ├── orm_analysis.py
│   ├── pathway_analysis.py
│   └── epigenetic_prediction.py
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
```

## Analysis Pipeline

### Step 1: Data Download & Preprocessing
```bash
python scripts/download_data.py
python scripts/data_preprocessing.py
```

### Step 2: ORM1/ORM2 Expression Analysis
```bash
python scripts/orm_analysis.py
```

### Step 3: Pathway Analysis
```bash
python scripts/pathway_analysis.py
```

### Step 4: Epigenetic Analysis
```bash
python scripts/epigenetic_prediction.py
```

### Step 5: Jupyter Notebooks
```bash
jupyter notebook
```
Then navigate to the `notebooks/` directory and run notebooks sequentially from 01 to 06.

## Expected Outputs

### Figures
- ORM1/ORM2 expression distributions (box plots, violin plots)
- Heatmaps of co-expressed genes
- Volcano plots for DEA
- Network diagrams
- UMAP/tSNE plots (if cell-level data available)

### Tables
- Differential expression results
- GO enrichment results
- KEGG pathway enrichment
- TFBS predictions
- Gene correlation matrix

### Reports
- HTML report with interactive visualizations
- Summary statistics
- Biological interpretation

## Key Analyses

### 1. Expression Pattern Analysis
Compare ORM1/ORM2 expression between:
- Psoriasis patients vs healthy controls
- Lesional vs non-lesional skin
- Different disease severity levels

### 2. Pathway Analysis
Identify biological processes and pathways enriched in ORM1/ORM2 co-expressed genes related to:
- Immune response
- Inflammatory signaling
- Acute phase response

### 3. Epigenetic Perspective
- Predict transcription factor binding sites (TFBS) in ORM1/ORM2 promoters
- Identify inflammatory TFs (NF-κB, STAT3, etc.)
- Infer chromatin accessibility changes

## Biological Context

**Orosomucoid (ORM)** also known as **alpha-1 acid glycoprotein (AGP)**:
- Acute phase protein produced primarily by hepatocytes
- Elevated in inflammatory conditions
- May have immunomodulatory functions
- Potential biomarker for psoriasis severity

## Author
Takahiro Yamamoto

## License
MIT License

## Citation
If you use this analysis pipeline, please cite:
- Federico et al., "Manually curated and harmonised transcriptomics datasets of psoriasis and atopic dermatitis." Sci Data 7, 153 (2020)
- Zenodo Dataset: https://zenodo.org/record/4009497

## References
1. Federico, A., et al. (2020). Manually curated and harmonised transcriptomics datasets of psoriasis and atopic dermatitis. Scientific Data, 7(1), 153.
2. Human Cell Atlas: https://data.humancellatlas.org/
3. ISD Atlas: https://isd-atlas.derma.meduniwien.ac.at/

## Contact & Issues
For questions or issues, please open a GitHub issue or contact the repository maintainer.