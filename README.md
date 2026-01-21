# Genomic Disease Association Analysis 🧬📊

### 🎓 MSc in Bioinformatics Project
This repository contains a statistical pipeline developed for the **"Programming in Python"** course as part of my **Master’s in Bioinformatics**. The project focuses on identifying significant associations between genomic data and disease categories.

### 🔍 Project Overview
The tool integrates **FlyBase** genomic annotations with the **Mondo Disease Ontology** to perform an enrichment analysis using structured biological data.

- **Ontology Integration**: Automatically maps Disease Ontology IDs (DOID) from FlyBase to higher-level Mondo categories using a tree-based traversal.
- **Statistical Significance**: Implements **Fisher's Exact Test** to evaluate the correlation between gene qualifiers and disease categories.
- **Advanced Statistics (Bonus Features)**: 
    - Calculates **Odds Ratio**, **Expected Values**, and **Fold Change** to measure the effect size of associations.
    - Implements **Multiple Testing Correction** using the **Benjamini-Hochberg (FDR)** method to ensure statistical reliability across thousands of comparisons.

### 🛠 Tech Stack
- **Python 3.x**
- **Scientific Libraries**: `SciPy` (statistical testing), `Statsmodels` (p-value correction).
- **Data Formats**: JSON (Mondo Ontology) and TSV (FlyBase Annotations).

### 🚀 How to Run
1. Ensure the following files are in the same directory:
   - `PROJECT_2_Genomic-Disease-Association-FisherTest.py`
   - `mondo_utils.py`
   - `mondo.json`
   - `disease_model_annotations_fb_2025_02.tsv`
2. Execute the script:
```bash
python PROJECT_2_Genomic-Disease-Association-FisherTest.py

```

📈 Example Output

```text

31522
unique genes = 5543

============================================================
TASK 6 10 Lowest P-values
============================================================
1. Mondo Category: neurodegenerative disease
   DO Qualifier:   ameliorates
   P-value:        9.1381e-169
   Odds-ratio:     7.9560
   Observed A:     689
   Expected A:     312.3242
   Fold Change:    2.2060

```
