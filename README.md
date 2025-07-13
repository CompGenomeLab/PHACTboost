# PHACTboost

PHACTboost is a gradiatent boosting tree based classifier that combines PHACT scores with information from multiple sequence alignment, phylogenetic trees, and ancestral reconstruction. 

The results of comprehensive experiments on carefully constructed sets of variants demonstrated that PHACTboost can outperform 40 prevalent pathogenicity predictors reported in the dbNSFP, including conventional tools, meta-predictors, and deep learning-based approaches as well as state-of-the-art tools, AlphaMissense, EVE, and CPT-1. The superiority of PHACTboost over these methods was particularly evident in case of hard variants for which different pathogenicity predictors offered conflicting results. We provide predictions of 219 million missense variants over 20,191 proteins. PHACTboost can improve our understanding of genetic diseases and facilitate more accurate diagnoses.

# Input Data

PHACTboost uses [PHACT](https://github.com/CompGenomeLab/PHACT) scores and their different versions (or related components) as features. The other input feature groups consist of gene and sequence position-specific features from the phylogenetic tree, ancestral probability distributions to integrate the structural properties of the phylogenetic tree and MSA-based frequency calculations as input features. Additionally, PHACTboost uses amino acid classes to utilize amino acid properties.

# PHACTboost Feature Construction Pipeline

## Requirements

### R Packages
```r
library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(Peptides)
library(Biostrings)
```

### Input Files
- Phylogenetic tree file (`.nwk` format)
- Ancestral probabilities file (`.state` format)
- Multiple sequence alignment (FASTA format)
- Log file from phylogenetic analysis
- IQTREE output file (`.iqtree`)
- Amino acid scales file (`Data/aa_scales.RData`)

## Pipeline Steps

### Step 1: MSA Masking (`MSA_Masking.R`)

**Usage**:
```bash
Rscript MSA_Masking.R <input_fasta> <uniprot_id> <output_folder>
```

**Parameters**:
- `input_fasta`: Raw FASTA file with protein sequences
- `uniprot_id`: UniProt identifier for the protein
- `output_folder`: Directory to save masked MSA

**Output**: Masked MSA file (`*_masked_msa.fasta`)


### Step 2: Feature Construction (`Main.R`)


**Usage**:
```bash
Rscript Main.R <tree_file> <ancestral_probs_file> <masked_msa_file> <log_file> <iqtree_file> <uniprot_id> <human_id> <parameters> <save_path> <aa_scales_file>
```

**Parameters**:
- `tree_file`: Phylogenetic tree file (`.nwk`)
- `ancestral_probs_file`: Ancestral probabilities file (`.state`)
- `masked_msa_file`: Masked MSA file from Step 1
- `log_file`: Log file from phylogenetic analysis (`.log`)
- `iqtree_file`: IQTREE output file (`.iqtree`)
- `uniprot_id`: UniProt identifier
- `human_id`: Human sequence identifier
- `parameters`: PHACT parameters (e.g., "CountNodes_3")
- `save_path`: Directory to save output files
- `aa_scales_file`: Amino acid scales file

**Output Files**:
- `*_scores.RData`: PHACT scores for each position
- `*_ml_features.RData`: Machine learning features
- `*_protscale_scores.RData`: Protein scale scores
- `*_protein_level_features.RData`: Protein-level statistics


### Step 3: Get Input Features (`get_input_features.R`)

**Usage**:
```bash
Rscript get_input_features.R <uniprot_id> <masked_msa_file>
```

**Parameters**:
- `uniprot_id`: UniProt identifier
- `masked_msa_file`: Masked MSA file from Step 1

**Output**: `input_features/<uniprot_id>.RData` - Final feature matrix

# Variant Set Construction

ClinVar (main data): variant_summary.txt dated 23.02.2023. (the closest version: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2023-02.txt.gz)

gnomAD: https://gnomad.broadinstitute.org/downloads#v3

HSV: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/homo_sapiens_variation.txt.gz

List of all variants added to ClinVar from 2015-2019 to construct training set:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2019 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2018 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2017 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2016 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2015 (11 files)

The GRCh38 coordinates reported in ClinVar, gnomAD and HSV are mapped to protein positions by using the pipeline reported in Map2Prot folder.

# Citing this work
Dereli, O.\*, Kuru, N.\*, Akkoyun, E., Bircan, A., Tastan, O.#, & Adebali, O.# (2024).  
*Molecular Biology and Evolution*, msae136.  
\* Co-first authors  
\# Co-corresponding authors 

https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msae136/7700170


# Data availability
All data obtained during this study is shared as Supplementary Material. 
Please see the [PHACTboost_manuscript](https://github.com/CompGenomeLab/PHACTboost_manuscript) folder for training and test sets of PHACTboost, as well as to reproduce the figures in the PHACTboost manuscript.

The entire prediction scores for PHACTboost and PHACT are also made available at [Zenodo](https://doi.org/10.5281/zenodo.11281312).

# Acknowledgement
This work was supported by the Health Institutes of Turkey (TUSEB) (Project no: 4587 to O.A.) and EMBO Installation Grant (Project no: 4163 to O.A.) funded by the Scientific and Technological Research Council of Turkey (TÜBİTAK). Most of the numerical calculations reported in this paper were performed at the High Performance and Grid Computing Center (TRUBA resources) of TÜBİTAK. The TOSUN cluster at Sabanci University was also used for computational analyses. We also want to thank Nehircan Özdemir for his art illustration of the PHACTboost approach.
