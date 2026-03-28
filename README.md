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

### UniProt
- Source: https://www.uniprot.org/
- Used to download pathogen proteomes
- Assembly ID required for matching

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
```

Activate the environment:

**Windows**
```bash
.\.venv\Scripts\Activate
```

**Linux / Mac**
```bash
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage and Description

Run scripts in this order:

```bash
python scripts/python/1_IEDB_Wrangling.py
python scripts/python/1_Pathogen_ref_wrangling.py
python scripts/python/2_Pathogen_fasta_retrivial.py
python scripts/python/3_Protein_meta_data.py
python scripts/python/4_Perfect_match_2_0.py
```

---

## 1. IEDB Wrangling (`1_IEDB_Wrangling.py`)

**Description**  
Filters epitopes by length (12–25 amino acids) and removes those with modified residues.  
To reduce redundancy, nested epitopes are identified by generating all possible 9-mers and removing shorter epitopes contained within longer ones.

**Input**
- Autoimmune epitopes from IEDB (TSV file)

**Output**
- Table with unique autoimmune epitopes  
- File: `data/intermediate/wrangled_IEDB.csv`

---

## 2. Pathogen Reference Wrangling (`1_Pathogen_ref_wrangling.py`)

**Description**  
Cleans and filters pathogen reference data. UniProt proteome data is cross-referenced with NCBI pathogen isolates using assembly IDs to identify human pathogenic organisms.

**Input**
- NCBI pathogen isolate data  
- UniProt proteome data  

**Output**
- List of human-pathogenic proteome IDs  
- File: `data/intermediate/proteome_ids.txt`

---

## 3. Proteome Retrieval (`2_Pathogen_fasta_retrivial.py`)

**Description**  
Uses the UniProt REST API to retrieve full FASTA proteomes for selected bacteria.

**Input**
- Proteome ID list from previous step  

**Output**
- FASTA files containing bacterial proteomes  
- File: `data/raw/all_fastas.fasta`

---

## 4. Protein Metadata Extraction (`3_Protein_meta_data.py`)

**Description**  
Extracts metadata from FASTA files, including organism, strain, annotation, gene name, sequence, and protein ID.

**Input**
- FASTA proteomes  

**Output**
- Table of protein metadata  
- File: `data/intermediate/wrangled_all_pathogen_prots.csv`

---

## 5. Epitope Matching (`4_Perfect_match_2_0.py`)

**Description**  
Performs epitope-protein matching using a sliding window approach (≥9-mers).  
Matching is performed using the Aho-Corasick algorithm, enabling efficient multi-pattern searching.

Only the longest match per epitope-pathogen-protein pair is retained. Unmatched proteins are also included and labeled.

**Input**
- Pathogen protein table  
- Autoimmune epitope table  

**Output**
- Table of matches and non-matches with metadata  
- File: `data/processed/perfect_matches_2_0.csv`

---
