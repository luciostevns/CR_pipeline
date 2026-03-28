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

---

## Usage and Description

Run scripts in this order:
- python 1_IEDB_Wrangling.py
- python 1_Pathogen_ref_wrangling.py
- python 2_Pathogen_fasta_retrivial.py
- python 3_Protein_meta_data.py
- python 4_Perfect_match_2_0.py

# 1_IEDB_Wrangling.py

Description:
This will filter epitope of being in specific length range (12-25) and having modified residues. to reduce redundancy in epitope list the script will check if smaller epitopes exist within other longer epitopes (nested epitopes), it will do this by creating all possible 9mers for each epitope if 2 epitopes are identical it will remove the parent epitope of the 9mer that is smallest. It will then write a .csv file, with chosen columns as columns and each row being a 9mer.

Input:
- Autoimmune epitopes sourced from IEDB.org as .tsv file

Output:
- table with unique autoimmune epitopes
- output file: "../Data/wrangled_IEDB.csv"

# 1_Pathogen_ref_wrangling

Description:
It will clean and filter data. It will take the uniport bacterial bulk data and check how many of them have been confirmed to be pathogenic against humans, it will do this with the ncbi isolates data by cross-referencing their assembly ID.

Input:
- ncbi.nlm.nih.gov/pathogens/isolates for confirmed pathogenicity and uniport.org data containing proteomes of choice.

Output:
- List of bacterial proteome IDs that is human pathogenic.
output file: "../Data/proteome_ids.txt"

# 2_Pathogen_fasta_retrivial

Description:
A rest api to retrive the full fasta proteomes of the bacteria.

Input:
- Output from 1_Pathogen_ref_wrangling (List of proteome IDs)

Output:
- 1 fasta file with all bacterial proteomes
- output file: "all_fastas.fasta"

# 3_Protein_meta_data

Description:
This will in essence just get all the meta data from the bacteria included. It will be each protein in the fastafiles metadata, so which Genus-species-strain each protein comes from, protein annotation, gene name, sequence and protein_ID.

Input:
- Output from 2_Pathogen_fasta_retrivial (bacterial proteomes fasta)

Output:
- Table with protein as rows and metadata as columns
- output file: "wrangled_all_pathogen_prots.csv"

# 4_Perfect_match_2_0 (Weird script name i know)

Description:
This is script does the main purpose of the project. The matching utilizes to main methods, to say it shortly it looks to see if the epitopes matches 1 to 1 with the pathogenic proteins. But since my epitopes can be up to 25 amino acids long, i utilize a sliding window approach that creates all possible sub epitopes down to 9mers, this results in many sub epitopes needed to be matched to a lot of proteins, so to do this within reasonable time, it utilizes an Aho-Corasick algorithm, this is a bit complicated but it sorta created a tree-structure where each branch is a sub epitope and if two sub-epitopes share a prefix they will share a branch until the differ where they branch off. there is more to it if you want to know more look at internet. Then the matching takes place where all sub-epitopes is compared to all protein, creating a lot of possible redundancy. We only want the longest match, meaning we only want 1 match per epitope, since one epitope might get matches on multiple of it sub epitopes, so we only keep the longest match per epitope-pathogen-protein pair. We also keep the pathogen proteins that did not match but flag them.

Input:
- output from 3_Protein_meta_data (pathogen protein table) and output from 1_IEDB_Wrangling.py (Autoimmune epitope table)

Output:
- Table where the rows is of each pathogen-protein-epitope match or pathogen-protein unmatched (so pathogen proteins not matched to autoimmune epitope is kept but labelled as such), with all metadata from the initial input tables.
- output file: "perfect_matches_2_0.csv"
