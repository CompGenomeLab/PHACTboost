load("/Users/nurdankuru/Desktop/PHACTboost/Train_gnomAD_Shared.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/Test_gnomAD_Shared.RData")

all <- rbind(train, test)

load("/Users/nurdankuru/Desktop/PHACTboost/dbNSFP_44a_TestSet.RData")
dbnsfp <- dbnsfp_testset

s <- colnames(dbnsfp)
tools <- colnames(dbnsfp)[grep("_rankscore", colnames(dbnsfp))]
check <- dbnsfp[,match(tools, colnames(dbnsfp))]

tools_considering_AF <- c("MetaSVM_rankscore", "MetaLR_rankscore", "MetaRNN_rankscore",
                          "BayesDel_addAF_rankscore", "BayesDel_noAF_rankscore", "M-CAP_rankscore",
                          "MutPred_rankscore", "ClinPred_rankscore")

check <- check[-match(tools_considering_AF, colnames(check))]

dis <- c()
keep <- c()
for (i in 1:length(check$SIFT_converted_rankscore)){
  print(i)
  vals <- as.numeric(check[i,which(check[i,]!=".")])
  
  if (length(vals)>=8){
    t <- length(vals)
    thr <- ceiling(t/2)-2
    l1 <- which(vals>0.5)
    l2 <- which(vals<0.5)
    if (length(l1)>=thr && length(l2)>=thr){
      keep <- c(keep, i)
      dis <- rbind(dis, c(length(l1), length(l2)))
    }
  }
}

vars <- dbnsfp$prot_vars[keep]

test <- test[match(vars, test$prot_vars),]
dbnsfp <- dbnsfp[match(vars, dbnsfp$prot_vars),]
dbnsfp_sub <- check[keep,]

print(sprintf("The variants of test set and dbNSFP matches? %s", length(sum(test$prot_vars==dbnsfp$prot_vars))))

hard_cases <- cbind(test$prot_vars, test$variant_info, test$`prediction$test_prediction`,
                    dbnsfp_sub)
colnames(hard_cases)[1:3] <- c("prot_vars", "variant_info", "PHACTboost")

print(table(hard_cases$variant_info))

save(hard_cases, file = "Hard_Cases.RData")
