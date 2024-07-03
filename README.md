# PHACTboost

PHACTboost is a gradiatent boosting tree based classifier that combines PHACT scores with information from multiple sequence alignment, phylogenetic trees, and ancestral reconstruction. 

The results of comprehensive experiments on carefully constructed sets of variants demonstrated that PHACTboost can outperform 40 prevalent pathogenicity predictors reported in the dbNSFP, including conventional tools, meta-predictors, and deep learning-based approaches as well as state-of-the-art tools, AlphaMissense, EVE, and CPT-1. The superiority of PHACTboost over these methods was particularly evident in case of hard variants for which different pathogenicity predictors offered conflicting results. We provide predictions of 219 million missense variants over 20,191 proteins. PHACTboost can improve our understanding of genetic diseases and facilitate more accurate diagnoses.

# Input Data

PHACTboost uses [PHACT](https://github.com/CompGenomeLab/PHACT) scores and their different versions (or related components) as features. The other input feature groups consist of gene and sequence position-specific features from the phylogenetic tree, ancestral probability distributions to integrate the structural properties of the phylogenetic tree and MSA-based frequency calculations as input features. Additionally, PHACTboost uses amino acid classes to utilize amino acid properties.


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
