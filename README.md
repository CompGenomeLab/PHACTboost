# PHACTboost

PHACTboost is a gradiatent boosting tree based classifier that combines PHACT scores with information from multiple sequence alignment, phylogenetic trees, and ancestral reconstruction. 

The results of comprehensive experiments on carefully constructed sets of variants demonstrated that PHACTboost can outperform 40 prevalent pathogenicity predictors reported in the dbNSFP, including conventional tools, meta-predictors, and deep learning-based approaches as well as state-of-the-art tools, AlphaMissense, EVE, and CPT-1. The superiority of PHACTboost over these methods was particularly evident in case of hard variants for which different pathogenicity predictors offered conflicting results. We provide predictions of 219 million missense variants over 20,191 proteins. PHACTboost can improve our understanding of genetic diseases and facilitate more accurate diagnoses.

## Quick Start

PHACTboost offers two usage modes:

### Option 1: Use Pre-trained Model (Recommended)
For making predictions on new variants using our pre-trained model:

**Prerequisites:**
1. Download the training data from [Google Drive](https://drive.google.com/file/d/1e433sHlmEzwSes858FngJbo3kP-5xg18/view?usp=drive_link)
2. Place the downloaded training data in the `Data/` directory

```bash
# Go to prediction directory
cd PHACTboost_Prediction/
Rscript PHACTboost_Prediction.R <ids> <codon_info_path> <train_path> <input_features_path> FinalModel/PHACTboost_model.txt
```

**Parameters:**
- `ids`: Vector of UniProt IDs to predict
- `codon_info_path`: Path to codon information file (`.xlsx`) (default: `Data/Codon.xlsx`)
- `train_path`: Path to training data for feature scaling reference (default: `Data/TrainingSet.RData`)
- `input_features_path`: Path to directory containing input features (`.RData` files) for the specific positions/variants you want to predict
- `final_model_path`: Path to the trained LightGBM model (default: `FinalModel/PHACTboost_model.txt`)

### Option 2: Train Your Own Model
For training PHACTboost from scratch:
```bash
# Go to machine learning directory
cd MachineLearning/
sbatch bash.sh
# or run directly
Rscript lgb.R <replication> <parameter_choice> <result_path>
```

**Parameters:**
- `replication`: Replication number for reproducibility
- `parameter_choice`: PHACT parameter choice (e.g., "CountNodes_3")
- `result_path`: Directory to save results

## Input Data

PHACTboost uses [PHACT](https://github.com/CompGenomeLab/PHACT) scores and their different versions (or related components) as features. The other input feature groups consist of gene and sequence position-specific features from the phylogenetic tree, ancestral probability distributions to integrate the structural properties of the phylogenetic tree and MSA-based frequency calculations as input features. Additionally, PHACTboost uses amino acid classes to utilize amino acid properties.

## PHACTboost Feature Construction Pipeline

This section describes how to construct PHACT features from scratch. This is required for both training new models and making predictions.

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
- Multiple sequence alignment (FASTA format)
- Amino acid scales file (`Data/aa_scales.RData`)

### Training Data
Due to size constraints, the training data is hosted externally:
- **Training Data**: [Download from Google Drive](https://drive.google.com/file/d/1e433sHlmEzwSes858FngJbo3kP-5xg18/view?usp=drive_link)
- Place the downloaded training data in the `Data/` directory

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

### Step 1.5: Re-run Ancestral State Reconstruction (ASR)

**Important**: After masking the MSA, you need to re-run the ancestral state reconstruction because the masked alignment changes the phylogenetic tree structure.

**Usage**:
```bash
iqtree2 -s ${masked_fasta} -te ${file_nwk} -m Data/vals.txt+R4 -asr --prefix ${id}_masked --safe
```

**Parameters**:
- `${masked_fasta}`: Masked MSA file from Step 1
- `${file_nwk}`: Original phylogenetic tree file (`.nwk`)
- `${id}`: UniProt identifier
- `Data/vals.txt`: Amino acid substitution model for ASR

**Output**: 
- `${id}_masked.state`: Ancestral probabilities file (needed for Step 2)
- `${id}_masked.treefile`: Updated phylogenetic tree
- `${id}_masked.iqtree`: IQTREE output file (needed for Step 2)
- `${id}_masked.log`: Log file (needed for Step 2)


### Step 2: Feature Construction (`Main.R`)


**Usage**:
```bash
Rscript Main.R <tree_file> <ancestral_probs_file> <masked_msa_file> <log_file> <iqtree_file> <uniprot_id> <human_id> <parameters> <save_path> <aa_scales_file>
```

**Parameters**:
- `tree_file`: Updated phylogenetic tree file from ASR analysis in Step 1.5 (`.treefile`)
- `ancestral_probs_file`: Ancestral probabilities file from ASR analysis in Step 1.5 (`.state`)
- `masked_msa_file`: Masked MSA file from Step 1
- `log_file`: Log file from ASR analysis in Step 1.5 (`.log`)
- `iqtree_file`: IQTREE output file from ASR analysis in Step 1.5 (`.iqtree`)
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
Rscript get_input_features.R <uniprot_id> <masked_msa_file> <data_path> <save_path> [aa_scales_path]
```

**Parameters**:
- `uniprot_id`: UniProt identifier
- `masked_msa_file`: Masked MSA file from Step 1
- `data_path`: Path to directory containing ML features (from Step 2)
- `save_path`: Directory to save output features (default: `input_features`)
- `aa_scales_path`: Path to amino acid scales file (default: `Data/aa_scales.RData`)

**Output**: `<save_path>/<uniprot_id>.RData` - Final feature matrix

## Model Training (MachineLearning/)

This section describes how to train PHACTboost models from scratch.

## Requirements
```r
library(lightgbm)
library(AUC)
```

## Usage
```bash
# Submit SLURM job
sbatch bash.sh

# Or run directly
Rscript lgb.R <replication> <parameter_choice> <result_path>
```

## Output
- Trained LightGBM model (`.txt` format)
- Cross-validation results
- Best hyperparameters
- Performance metrics (train/test AUC)

## Prediction (PHACTboost_Prediction/)

This section describes how to make predictions using the pre-trained PHACTboost model.

## Requirements
```r
library(lightgbm)
library(AUC)
library(readxl)
```

## Usage
```bash
Rscript PHACTboost_Prediction.R <ids> <codon_info_path> <train_path> <input_features_path> FinalModel/PHACTboost_model.txt
```

## Parameters
- `ids`: Vector of UniProt IDs to predict
- `codon_info_path`: Path to codon information file (`.xlsx`) (default: `Data/Codon.xlsx`)
- `train_path`: Path to training data for feature scaling reference (default: `Data/TrainingSet.RData`)
- `input_features_path`: Path to directory containing input features (`.RData` files) for the specific positions/variants you want to predict
- `final_model_path`: Path to the trained LightGBM model (default: `FinalModel/PHACTboost_model.txt`)

## Output
- `PHACTboost_<id>.RData`: Prediction results with PHACTboost scores

## Variant Set Construction

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

## Citing this work
Dereli, O.\*, Kuru, N.\*, Akkoyun, E., Bircan, A., Tastan, O.#, & Adebali, O.# (2024).  
*Molecular Biology and Evolution*, msae136.  
\* These authors contributed equally to this work. Cofirst authors O.D. and N.K. have the right to list themselves first in the
author order on their CVs.
\# Co-corresponding authors 

https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msae136/7700170


## Data availability
All data obtained during this study is shared as Supplementary Material. 

### Training and Test Sets
Complete feature matrices for training and testing PHACTboost models:
- **Training Set**: [Download from Google Drive](https://drive.google.com/file/d/1e433sHlmEzwSes858FngJbo3kP-5xg18/view?usp=drive_link)
- **Test Set**: [Download from Google Drive](https://drive.google.com/file/d/1LfHYKRu3ro5HyaGUeLb5qh9ajtMhEzYS/view?usp=drive_link)

These datasets contain all computed features and are ready for model training and evaluation.

### Additional Resources
Please see the [PHACTboost_manuscript](https://github.com/CompGenomeLab/PHACTboost_manuscript) folder for training and test sets of PHACTboost, as well as to reproduce the figures in the PHACTboost manuscript.

The entire prediction scores for PHACTboost and PHACT are also made available at [Zenodo](https://doi.org/10.5281/zenodo.11281312). This dataset contains predictions for all human proteins at SNV (Single Nucleotide Variant) and amino acid substitution level.

## Acknowledgement
This work was supported by the Health Institutes of Turkey (TUSEB) (Project no: 4587 to O.A.) and EMBO Installation Grant (Project no: 4163 to O.A.) funded by the Scientific and Technological Research Council of Turkey (TÜBİTAK). Most of the numerical calculations reported in this paper were performed at the High Performance and Grid Computing Center (TRUBA resources) of TÜBİTAK. The TOSUN cluster at Sabanci University was also used for computational analyses. We also want to thank Nehircan Özdemir for his art illustration of the PHACTboost approach.
