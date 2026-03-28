# Computational Prediction of T-cell Cross-Reactivity

## Overview

This project aims to identify potential sequence similarity between known autoimmune T-cell epitopes and pathogen proteomes. The goal is to investigate whether molecular mimicry may contribute to autoimmune responses.

The pipeline integrates data from IEDB, UniProt, and NCBI to detect shared peptide motifs between human autoantigens and pathogenic proteins.

---

## Data Sources

### IEDB (Immune Epitope Database)
- Source: https://www.iedb.org/
- Contains experimentally validated T-cell epitopes
- Required fields:
  - Assay ID
  - Epitope sequence
  - Protein source
  - Protein ID
  - Disease
  - Disease stage
  - MHC restriction
  - Epitope start/end positions
  - Modified residues

---

### UniProt
- Source: https://www.uniprot.org/
- Used to download pathogen proteomes
- Assembly ID required for matching

---

### NCBI Pathogen Database
- Source: https://www.ncbi.nlm.nih.gov/pathogens/
- Used to verify pathogenicity of selected organisms

---

## Project Structure

```
CR_pipeline/
├── data/
│   ├── raw/
│   ├── intermediate/
│   └── processed/
├── results/
│   ├── figures/
│   └── tables/
├── scripts/
│   ├── python/
│   └── r/
├── requirements.txt
└── README.md
```

---

## Setup

Clone the repository and create a virtual environment:

```bash
python -m venv .venv
