load("/Users/nurdankuru/Desktop/dbNSFP_43a_AllVariantSet_AllColumns.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/Train_gnomAD_Shared.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/Test_gnomAD_Shared.RData")

all <- rbind(train, test)
rm(train, test)

data$chr_var_with_aa <- paste(data$`#chr`, "-", data$`pos(1-based)`, data$ref, ">", data$alt, "-",
                              data$aaref, data$aaalt, sep = "")
all$chr_var_with_aa <- paste(all$chr_vars, "-", all$Ref_AA, all$Alt_AA, sep = "")

common_vars <- intersect(data$chr_var_with_aa, all$chr_var_with_aa)

dbnsfp_allvar <- data[match(common_vars, data$chr_var_with_aa),]
all_sub <- all[match(common_vars, all$chr_var_with_aa),]

print(sprintf("Number of common variants is %s = %s", sum(dbnsfp_allvar$chr_var_with_aa==all_sub$chr_var_with_aa),
              length(dbnsfp_allvar$chr_var_with_aa)))

dbnsfp_allvar$prot_vars <- all_sub$prot_vars
dbnsfp_allvar$variant_info <- all_sub$variant_info
dbnsfp_allvar$UNIPROTKB <- all_sub$UNIPROTKB

tools <- grep("rankscore", colnames(dbnsfp_allvar))
tools_considering_AF <- c("MetaSVM_rankscore", "MetaLR_rankscore", "MetaRNN_rankscore",
                          "BayesDel_addAF_rankscore", "BayesDel_noAF_rankscore", "M-CAP_rankscore",
                          "MutPred_rankscore", "ClinPred_rankscore", "Eigen-raw_coding_rankscore")
tools_lessnumvar <- c("LINSIGHT_rankscore")

tools <- tools[-match(tools_considering_AF, colnames(dbnsfp_allvar)[tools])]
tools <- tools[-match(tools_lessnumvar, colnames(dbnsfp_allvar)[tools])]

dbnsfp_allvar <- dbnsfp_allvar[, c(1:13, tools, 14:37, 786:789)]

save(dbnsfp_allvar, file="dbNSFP_44a_CompleteSet.RData")

load("Test_gnomAD_Shared.RData")

common_vars_test <- intersect(test$prot_vars, dbnsfp_allvar$prot_vars)
dbnsfp_testset <- dbnsfp_allvar[match(common_vars_test, dbnsfp_allvar$prot_vars),]

save(dbnsfp_testset, file="dbNSFP_44a_TestSet.RData")
