library(data.table)
library(stringr)
library(readxl)
library(dplyr)

clinvar_bef <- "variant_summary.txt"

clinvar_bef <- fread(file = clinvar_bef, sep = "\t", header = T, data.table = F)
clinvar_bef_38 <- clinvar_bef[clinvar_bef$Assembly == "GRCh38", ]
clinvar_bef_38 <- clinvar_bef_38[clinvar_bef_38$Type == "single nucleotide variant", ]
clinvar_bef_38$chr_vars <- paste0(clinvar_bef_38$Chromosome, "-", clinvar_bef_38$Start, clinvar_bef_38$ReferenceAlleleVCF, ">", clinvar_bef_38$AlternateAlleleVCF)

save(clinvar_bef_38, file="Clinvar_bef_38.RData")

#load("Clinvar_bef_38.RData")
load("ClinVar_Mapped.RData")

clinvar <- clinvar[which(clinvar$mutation_type=="non_synonymous"),]
chr_vars <- paste(unlist(clinvar[,1]), "-", unlist(clinvar[,2]), unlist(clinvar[,4]),
                  ">", unlist(clinvar[,5]), sep = "")
clinvar$chr_vars <- chr_vars
UNIPROTKB <- str_split(unlist(clinvar$original_protein_id), "\\|", simplify = T)
clinvar$UNIPROTKB <- unlist(UNIPROTKB[,2])

clinvar$prot_vars <- paste(unlist(clinvar$UNIPROTKB), "-", unlist(clinvar$alt_protein_position), 
                           unlist(clinvar$original_aa), "/", unlist(clinvar$alt_aa), sep = "")

print(length(intersect(clinvar$chr_vars, clinvar_bef_38$chr_vars)))

print("PART 1 - ClinVar Data")
save(clinvar, file="clinvar.RData")

print(clinvar[1:3,])

clinvar_nonsy <- merge(clinvar_bef_38, clinvar[,c("prot_vars", "chr_vars", "mutation_type", "UNIPROTKB")],
                       by = "chr_vars", all.x = T)
clinvar_nonsy <- clinvar_nonsy[-which(is.na(clinvar_nonsy$UNIPROTKB)),]
save(clinvar_nonsy, file="ClinVar_NonSyn.RData")

#########################################
##### Read gnomAD File #####

load("gnomad_v3_2.RData")
load("gnomAD_Mapped.RData")

print("PART 2 - gnomAD Data")

gnomad$chr_vars <- paste(unlist(gnomad[,1]), "-", unlist(gnomad[,2]), unlist(gnomad[,4]), 
                         ">", unlist(gnomad[,5]), sep = "")
gnomad <- gnomad[which(gnomad$mutation_type=="non_synonymous"),]

UNIPROTKB <- str_split(unlist(gnomad$original_protein_id), "\\|", simplify = T)
gnomad$UNIPROTKB <- unlist(UNIPROTKB[,2])

gnomad$prot_vars <- paste(unlist(gnomad$UNIPROTKB), "-", unlist(gnomad$alt_protein_position), 
                           unlist(gnomad$original_aa), "/", unlist(gnomad$alt_aa), sep = "")

gnomad_nonsy <- merge(gnomad_v3_2, gnomad[,c("prot_vars", "chr_vars", "mutation_type", "UNIPROTKB")],
                       by = "chr_vars", all.x = T)
gnomad_nonsy <- gnomad_nonsy[-which(is.na(gnomad_nonsy$UNIPROTKB)),]
save(gnomad_nonsy, file="gnomAD_Nonsyn.RData")

rm(clinvar, clinvar_bef_38, gnomad, gnomad_v3_2, UNIPROTKB)

data_all <- clinvar_nonsy
data_all <- merge(data_all, gnomad_nonsy[,c("AN", "AC", "AF", "chr_vars", "prot_vars")], by = "chr_vars", all.x = T)

not_exist_in_clinvar <- match(setdiff(gnomad_nonsy$chr_vars, data_all$chr_vars), gnomad_nonsy$chr_vars)
neutrals_from_gnomad <- which(as.numeric(gnomad_nonsy$AF)>=0.01)
neutrals <- gnomad_nonsy[intersect(not_exist_in_clinvar, neutrals_from_gnomad),]

print(dim(data_all))
print(dim(neutrals))

save(neutrals, file="NeutralOnlygnomAD.RData")
save(data_all, file="ClinVar_withAF.RData")

rm(gnomad_nonsy, clinvar_nonsy)

colnames(data_all)[36] <- "prot_vars"
colnames(data_all)[42] <- "prot_vars_gnomad"

################################################################## 
##### Get Submission Dates Data from "submisson_summary.txt" #####

load("AllInfo.RData")
vars <- data_all$VariationID
dates <- as.data.frame(dates)
sub_aft <- match(vars, dates$V1)

clinvar_with_dates <- cbind(data_all, dates$V2[sub_aft])
colnames(clinvar_with_dates)[43] <- "Dates"

###########################################
##### Assign Star Info & Variant Type #####

star_info <- matrix(-1, length(clinvar_with_dates$chr_vars),1)
star_info[which(clinvar_with_dates$ReviewStatus== "criteria provided, multiple submitters, no conflicts")] <- 2
star_info[which(clinvar_with_dates$ReviewStatus== "criteria provided, single submitter")] <- 1
star_info[which(clinvar_with_dates$ReviewStatus== "no assertion criteria provided")] <- 0
star_info[which(clinvar_with_dates$ReviewStatus== "no assertion provided")] <- 0
star_info[which(clinvar_with_dates$ReviewStatus== "no interpretation for the single variant")] <- 0
star_info[which(clinvar_with_dates$ReviewStatus== "practice guideline")] <- 4
star_info[which(clinvar_with_dates$ReviewStatus== "reviewed by expert panel")] <- 3

print(table(clinvar_with_dates$ReviewStatus))

clinvar_with_dates <- cbind(star_info, clinvar_with_dates)

save(clinvar_with_dates, file = "ClinVar_WithDates.RData")

print(dim(clinvar_with_dates))
print("PART 3 - Star Assigned")

##########################################################
##### Neutrals Existing at gnomAD but not in ClinVar #####

add_to_clinvar <- cbind(matrix(0, length(neutrals$CHROM), 39), neutrals[,c("AN", "AC", "AF")], "-", "-")
colnames(add_to_clinvar) <- colnames(clinvar_with_dates)
add_to_clinvar$chr_vars <- neutrals$chr_vars
add_to_clinvar$prot_vars <- neutrals$prot_vars
add_to_clinvar$prot_vars_gnomad <- neutrals$prot_vars
add_to_clinvar$star_info <- 1
add_to_clinvar$ReviewStatus <- "criteria provided, single submitter"
add_to_clinvar$ClinicalSignificance <- "Benign"
add_to_clinvar$Chromosome <- neutrals$CHROM
add_to_clinvar$Start <- as.numeric(neutrals$POS)+1
add_to_clinvar$Stop <- as.numeric(neutrals$POS)+1
add_to_clinvar$ReferenceAlleleVCF <- neutrals$REF
add_to_clinvar$AlternateAlleleVCF <- neutrals$ALT
add_to_clinvar$ReferenceAllele <- neutrals$REF
add_to_clinvar$AlternateAllele <- neutrals$ALT
add_to_clinvar$Origin <- "germline"
add_to_clinvar$OriginSimple <- "germline"
add_to_clinvar$Assembly <- "GRCh38"
add_to_clinvar$PhenotypeList <- "Benign - AF"
add_to_clinvar$PhenotypeIDS <- "-"
add_to_clinvar$RCVaccession <- "-"
add_to_clinvar$GeneID <- "-"
add_to_clinvar$GeneSymbol <- "-"
add_to_clinvar$VariationID <- "-"
add_to_clinvar$PositionVCF <- as.numeric(neutrals$POS)+1
add_to_clinvar$UNIPROTKB <- neutrals$UNIPROTKB
add_to_clinvar$mutation_type <- neutrals$mutation_type
add_to_clinvar$Dates <- "1900"
add_to_clinvar$Type <- "single nucleotide variant"
save(add_to_clinvar, file="add_to_clinvar.RData")

clinvar_with_dates <- rbind(clinvar_with_dates, add_to_clinvar)

print(dim(clinvar_with_dates))
save(clinvar_with_dates, file = "ClinVar_WithDates_and_gnomAD.RData")

print("PART 4 - Merged ClinVar&gnomAD")

clinvar_mapped <- clinvar_with_dates

##############################################
##### ELIM 1: Take P, LP, B, LB Variants #####

label1 <- "Pathogenic"
label2 <- "Benign"
label3 <- "Likely pathogenic"
label4 <- "Pathogenic/Likely pathogenic"
label5 <- "Benign/Likely benign"
label6 <- "Likely benign"

ll <- c(which(clinvar_mapped$ClinicalSignificance==label1), which(clinvar_mapped$ClinicalSignificance==label2),
        which(clinvar_mapped$ClinicalSignificance==label3), which(clinvar_mapped$ClinicalSignificance==label4),
        which(clinvar_mapped$ClinicalSignificance==label5), which(clinvar_mapped$ClinicalSignificance==label6))

clinvar_mapped <- clinvar_mapped[ll,]

save(clinvar_mapped, file = "ClinVar_Mapped_gnomad_elimVUS.RData")

print(dim(clinvar_mapped))
print("ELIM 1")

####################################################
##### ELIM 2: Eliminate Conf. Variants wrt HSV #####

load("missense_variation.RData")

mis <- merge(missense_variation, clinvar_mapped, by = "chr_vars", all.x = T)
mis <- mis[-which(is.na(mis$star_info)==1),]
save(mis, file="Subset_HSV_NewMapping_Merge.RData")
rm(missense_variation)

sub_missense <- mis
matched <- cbind(sub_missense$chr_vars, sub_missense[,7], sub_missense[,23])

uncertains <- matched[grep("ncertain", matched[,2]),]
matched <- matched[-grep("ncertain", matched[,2]),]
uncertains2 <- matched[intersect(grep("enign", matched[,2]), grep("athogen", matched[,3])),]
matched <- matched[-intersect(grep("enign", matched[,2]), grep("athogen", matched[,3])),]
uncertains3 <- matched[intersect(grep("enign", matched[,3]), grep("athogen", matched[,2])),]
matched <- matched[-intersect(grep("enign", matched[,3]), grep("athogen", matched[,2])),]
uncertains4 <- matched[grep("association", matched[,2]),]
matched <- matched[-grep("association", matched[,2]),]
uncertains5 <- matched[grep("onflicting", matched[,2]),]
matched <- matched[-grep("onflicting", matched[,2]),]
uncertains6 <- matched[grep("rotective", matched[,2]),]
matched <- matched[-grep("rotective", matched[,2]),]
uncertains7 <- matched[grep("rug resp", matched[,2]),]
matched <- matched[-grep("rug resp", matched[,2]),]
uncertains8 <- matched[grep("isk factor", matched[,2]),]
matched <- matched[-grep("isk factor", matched[,2]),]

problem_all <- rbind(uncertains, uncertains2, uncertains3, uncertains4, uncertains5, uncertains6, uncertains7,
                     uncertains8)
rm(uncertains, uncertains2, uncertains3, uncertains4, uncertains5, uncertains6, uncertains7, uncertains8)
prob <- intersect(grep("enign", matched[,2]), grep("athoge", matched[,3]))
prob <- intersect(grep("enign", matched[,3]), grep("athoge", matched[,2]))

print(dim(prob))

eliminate_because_of_HSV_vars <- unique(problem_all[,1])
dubles <- which(table(clinvar_mapped$chr_vars)>1)
for (i in names(dubles)){
  delete <- which(clinvar_mapped$chr_vars==i)
  keep_del <- clinvar_mapped[delete, ]
  
  clinvar_mapped <- clinvar_mapped[-delete,]
  
  keep_add <- matrix(0, 1, length(clinvar_mapped[1,]))
  for (i in 1:length(clinvar_mapped[1,])){
    keep_add[i] <- paste(unique(keep_del[,i]), sep = ";")
  }
  colnames(keep_add) <- colnames(clinvar_mapped)
  clinvar_mapped <- rbind(clinvar_mapped, keep_add)
}

eliminate_because_of_HSV <- clinvar_mapped[match(eliminate_because_of_HSV_vars, clinvar_mapped$chr_vars),]

clinvar_mapped <- clinvar_mapped[-match(eliminate_because_of_HSV_vars, clinvar_mapped$chr_vars),]

print(dim(clinvar_mapped))
save(clinvar_mapped, file = "ClinVar_Mapped_HSV_Eliminated.RData")

print("ELIM 2")

########################################
##### ELIM 3: Assign Variant Info #####

variant_info <- matrix(0, length(clinvar_mapped$star_info), 1)
variant_info[grep("enign", clinvar_mapped$ClinicalSignificance)] <- -1
variant_info[grep("athogen", clinvar_mapped$ClinicalSignificance)] <- 1
clinvar_mapped <- cbind(variant_info, clinvar_mapped)

########################################
##### ELIM 3: Eliminate Duplicates #####

clinvar_mapped_no_dup <- matrix(0, length(unique(clinvar_mapped$prot_vars)), length(clinvar_mapped[1,]))
k <- 1
ch <- 0
unq_vars <- unique(clinvar_mapped$prot_vars)
for (i in unq_vars){
  w <- which(clinvar_mapped$prot_vars==i)
  if (length(w)>1){
    line <- matrix(0,1,length(clinvar_mapped[1,]))
    for (ind in 1:length(clinvar_mapped[1,])){
      val <- unique(clinvar_mapped[w,ind])
      line[ind] <- paste(val, collapse = ";")
    }
  } else {
    line <- clinvar_mapped[w,]
  }
  
  clinvar_mapped_no_dup[k, ] <- unlist(line)
  k <- k + 1
}

colnames(clinvar_mapped_no_dup) <- colnames(clinvar_mapped)
clinvar_mapped_no_dup <- as.data.frame(clinvar_mapped_no_dup)

print(dim(clinvar_mapped_no_dup))
save(clinvar_mapped_no_dup, file="clinvar_mapped_no_dup_Final.RData")

print("ELIM 3")

############################################################
##### ELIM 4: Eliminate Incompatible Labelled Variants #####

el_incomp <- grep(";", clinvar_mapped_no_dup$variant_info)
xx <- clinvar_mapped_no_dup[el_incomp,]
if (length(el_incomp)>1){
  clinvar_mapped_no_dup <- clinvar_mapped_no_dup[-el_incomp,]
}
clinvar_mapped <- clinvar_mapped_no_dup

print(dim(clinvar_mapped))
save(clinvar_mapped, file = "ClinVar_Mapped_HSV_Eliminated_NoDup.RData")

print("ELIM 4")

#################################################
##### Determine Neutrals & Pathogens wrt AF #####

pathogens <- clinvar_mapped[grep("athogen", clinvar_mapped$ClinicalSignificance),]
neutrals <- clinvar_mapped[grep("enign", clinvar_mapped$ClinicalSignificance),]

nv <- grep(";", pathogens$AF)
AF_patogen <- pathogens$AF
for (i in nv){
  af <- unlist(strsplit(AF_patogen[i], ";"))
  if (length(which(af!="NA"))>=1){
    AF_patogen[i] <- max(as.numeric(af[which(af!="NA")]))
  } else {
    AF_patogen[i] <- "NA"
  }
}
AF_patogen[c(which(is.na(AF_patogen)), which(AF_patogen=="NA"))] <- 0
el_pat <- which(as.numeric(AF_patogen)>=0.005) #-23
pathogens <- pathogens[-el_pat,]

nv <- grep(";", neutrals$AF)
AF_neutral <- neutrals$AF
for (i in nv){
  af <- unlist(strsplit(AF_neutral[i], ";"))
  if (length(which(af!="NA"))>=1){
    AF_neutral[i] <- min(as.numeric(af[which(af!="NA")]))
  } else {
    AF_neutral[i] <- "NA"
  }
}
el1 <- which(is.na(AF_neutral))
el2 <- which(AF_neutral=="NA")
AF_neutral[unique(c(el1, el2))] <- 0
el_net <- which(as.numeric(AF_neutral)<0.01) #57388
el_net_all <- unique(c(el1, el2, el_net))
neutrals <- neutrals[-el_net_all,]

clinvar_mapped_final <- rbind(pathogens, neutrals)

print(dim(clinvar_mapped_final))
save(clinvar_mapped_final, file="ClinVar_Final.RData")

print("Determine Neutrals & Pathogens")

#####################################
##### TRAIN & TEST SET - Step 1 #####

load("All_Chr_Vars_2016_2020.RData")
colnames(chr_vars_all) <- c("Year", "ChrVars", "Type", "Stat", "Date")
chr_vars_all <- as.data.frame(chr_vars_all)
chr_vars_16_20 <- chr_vars_all
chr_vars_2020 <- chr_vars_16_20[which(chr_vars_16_20[,1]=="2020"),]
chr_vars_2019 <- chr_vars_16_20[which(chr_vars_16_20[,1]=="2019"),]
chr_vars_2018 <- chr_vars_16_20[which(chr_vars_16_20[,1]=="2018"),]
chr_vars_2017 <- chr_vars_16_20[which(chr_vars_16_20[,1]=="2017"),]
chr_vars_2016 <- chr_vars_16_20[which(chr_vars_16_20[,1]=="2016"),]
chr_vars_2015 <- chr_vars_16_20[which(chr_vars_16_20[,1]=="2015"),]
load("All_Chr_Vars_2015.RData")
chr_vars_2015 <- chr_vars_all
current <- clinvar_mapped_final$chr_vars


all_15_19 <- unique(c(unlist(chr_vars_2015[,2]), unlist(chr_vars_2016[,2]), unlist(chr_vars_2017[,2]), 
                      unlist(chr_vars_2018[,2]), unlist(chr_vars_2019[,2])))
all_15_19 <- unique(all_15_19)

ones <- clinvar_mapped_final[-grep(";", clinvar_mapped_final$chr_vars),]
twices <- clinvar_mapped_final[grep(";", clinvar_mapped_final$chr_vars),]

### PART1
before_19 <- intersect(all_15_19, ones$chr_vars)
train_part1 <- clinvar_mapped_final[match(before_19, clinvar_mapped_final$chr_vars),]

### PART 2
tr_ind <- c()
for (i in twices$chr_vars){
  i2 <- unlist(strsplit(i, ";"))
  if (sum(is.element(i2, all_15_19))>0){
    tr_ind <- c(tr_ind, i)
  }
}
train_part2 <- clinvar_mapped_final[match(tr_ind, clinvar_mapped_final$chr_vars),]

train_set <- rbind(train_part1, train_part2)
remainings <- setdiff(clinvar_mapped_final$chr_vars, train_set$chr_vars)
test_set <- clinvar_mapped_final[match(remainings, clinvar_mapped_final$chr_vars),]

print("Train & Test - Part 1")

#####################################
##### TRAIN & TEST SET - Step 2 #####

el_date_tr <- c()
el_date_ts <- c()
for (i in 1900:2019){
  el <- grep(i, train_set$Dates)
  el_date_tr <- c(el_date_tr, el)
  el2 <- grep(i, test_set$Dates)
  el_date_ts <- c(el_date_ts, el2)
}
el_date_tr <- unique(el_date_tr)
el_date_ts <- unique(el_date_ts)

train_final <- rbind(train_set, test_set[el_date_ts,])
test_final <- test_set[-el_date_ts,]

test_final <- test_final[-which(test_final$star_info==0),]
save(train_final, file="Train_2019_gnomAD_Train.RData")
save(test_final, file="Test_2019_gnomAD_Test.RData")

print("Train & Test - Part 2")


k <- c()
for (i in train_final$chr_vars){
  i <- unlist(strsplit(i, ";"))
  
  k <- c(k, sum(is.element(i, all_15_19)))
}

print(c("Number of common vars in 15-20 set ", sum(k)))

k <- c()
l <- 0
for (i in test_final$chr_vars){
  print(l)
  l <- l+1
  i <- unlist(strsplit(i, ";"))
  
  k <- c(k, sum(is.element(i, all_15_19)))
}

print(c("Number of common variants in 15-20 and test ", sum(k)))

