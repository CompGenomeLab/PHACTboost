# PHACTboost

PHACTboost is a gradiatent boosting tree based classifier that combines PHACT scores with information from multiple sequence alignment, phylogenetic trees, and ancestral reconstruction. Our comprehensive experiments on 20,191 human proteins encompassing 219 million amino acid substitution  demonstrate that PHACTboost overpeforms 40 available pathogenicity predictors  including the recent ones such as EVE and CPT-1. 

# Variant Set Construction

ClinVar (main data): https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2023-02.txt.gz
gnomAD: 
HSV: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/homo_sapiens_variation.txt.gz

List of all variants added to ClinVar from 2015-2019 to construct training set:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2019 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2018 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2017 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2016 (12 files)
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2015 (11 files)


# Citing this work

# Data availability
All data obtained during this study is shared on the following url. 

* Please see Publications/... folder to reproduce the figures in PHACTboost manuscript.

# Acknowledgement
This work was supported by the Health Institutes of Turkey (TUSEB) (Project no: 4587 to O.A.) and EMBO Installation Grant (Project no: 4163 to O.A.) funded by the Scientific and Technological Research Council of Turkey (TÜBİTAK). Most of the numerical calculations reported in this paper were performed at the High Performance and Grid Computing Center (TRUBA resources) of TÜBITAK. The TOSUN cluster at Sabanci University was also used for computational analyses. We also want to thank Nehircan Özdemir for his art illustration of the PHACTboost approach.
