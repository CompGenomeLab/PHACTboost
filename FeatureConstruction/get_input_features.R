#!/usr/bin/env Rscript
library(bio3d)
library(stringr)
library(Biostrings)
library(Peptides)

args <- commandArgs(trailingOnly = TRUE)
uniprot_id <- args[1]

save_path <- "input_features_New"

if (dir.exists(sprintf("%s", save_path)) == FALSE) {
  dir.create(sprintf("%s", save_path)) 
}

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

keep_gaps <- c()
parameters <- c("CountNodes_3")
aa_sift <- "  A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z   *   -"
aa_phact <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")

ids <- uniprot_id
#ids <- list.files("MLFeats")

for (id in ids){
  print(id)
  
  data_pathh <- sprintf("MLFeats")
  load(sprintf("%s/%s/%s_scores.RData", data_pathh,  id, id))
  load(sprintf("%s/%s/%s_ml_features.RData", data_pathh,  id, id))
  load(sprintf("%s/%s/%s_protscale_scores.RData", data_pathh, id, id))
  load(sprintf("%s/%s/%s_protein_level_features.RData", data_pathh, id, id))
  
  file_fasta <- paste("/cta/groups/adebali/msa_masking_nurdan/Alignments/", id, "_MaskedMSA.fasta", sep = "")
  # Read fasta file, MSA
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  human_seq <- msa[grep(id, row.names(msa)),]
  if (length(which(human_seq=="-"))>=1){
    keep_gaps <- which(human_seq=="-")
  }
  
  positions <- c(1:dim(scores[[sprintf("%s_wl_param_CountNodes_3", id)]])[1])
  
  ref_aa <- sapply(positions, function(x){ml_features[[x]]$ref_aa})
  names(ref_aa) <- NULL
  
  positions_all <- c(sapply(positions, function(x){rep(positions[x], 20)}))
  
  ref_aa_all <- c(sapply(positions, function(x){rep(ref_aa[x], 20)}))
  alt_aa_all <- rep(aa_phact, length(positions))
  
  data <- cbind(id, positions_all, ref_aa_all, alt_aa_all)
  colnames(data) <- c("UNIPROTKB", "Positions", "Ref_AA", "Alt_AA")
  
  data <- as.data.frame(data)
  
  positions <- as.numeric(data$Positions)
  s1 <- sapply(parameters, function(p){
    scores_wp_wl <- scores[[sprintf("%s_wl_param_%s", id, p)]]
    scores_wp_wol <- scores[[sprintf("%s_wol_param_%s", id, p)]]
    scores_no_gap <- scores[[sprintf("%s_nogap_param_%s", id, p)]]
    
    ref_aa <- match(data$Ref_AA, colnames(scores_wp_wl))
    alt_aa <- match(data$Alt_AA, colnames(scores_wp_wl))
    
    c1 <- scores_wp_wl[cbind(positions, ref_aa)]
    c2 <- scores_wp_wl[cbind(positions, alt_aa)]
    c3 <- scores_wp_wol[cbind(positions, ref_aa)]
    c4 <- scores_wp_wol[cbind(positions, alt_aa)]
    
    c5 <- scores_no_gap[cbind(positions, ref_aa)]
    c6 <- scores_no_gap[cbind(positions, alt_aa)]
    
    data[[sprintf("Phact_wl_param_%s_ref", p)]] <<- c1
    data[[sprintf("Phact_wl_param_%s_alt", p)]] <<- c2
    data[[sprintf("Phact_wol_param_%s_ref", p)]] <<- c3
    data[[sprintf("Phact_wol_param_%s_alt", p)]] <<- c4
    
    data[[sprintf("Phact_wl_param_%s_ratio", p)]] <<- c1 / c2
    data[[sprintf("Phact_wol_param_%s_ratio", p)]] <<- c3 / c4
    
    diversity <- scores[[sprintf("%s_diversity_%s", id, p)]]
    #  diversity <- 1 - log(diversity + 1e-15)/log(1e-15)
    data[[sprintf("Phact_diversity_param_%s", p)]] <<- diversity[positions]
    
    data[[sprintf("Phact_nogap_param_%s_ref", p)]] <<- c5
    data[[sprintf("Phact_nogap_param_%s_alt", p)]] <<- c6
    
    score_max_wl <- apply(scores_wp_wl[,2:21], 1, max)
    score_max_wol <- apply(scores_wp_wol[,2:21], 1, max)
    
    score_min_wl <-  apply(scores_wp_wl[,2:21], 1, min)
    score_min_wol <- apply(scores_wp_wol[,2:21], 1, min)
    
    score_sd_wl <- apply(scores_wp_wl[,2:21], 1, sd)
    score_sd_wol <- apply(scores_wp_wol[,2:21], 1, sd)
    
    score_avg_wl <- apply(scores_wp_wl[,2:21], 1, mean)
    score_avg_wol <- apply(scores_wp_wol[,2:21], 1, mean)
    
    data[[sprintf("Phact_wl_param_%s_max", p)]] <<- score_max_wl[positions]
    data[[sprintf("Phact_wol_param_%s_max", p)]] <<- score_max_wol[positions]
    
    data[[sprintf("Phact_wl_param_%s_min", p)]] <<- score_min_wl[positions]
    data[[sprintf("Phact_wol_param_%s_min", p)]] <<- score_min_wol[positions]
    
    data[[sprintf("Phact_wl_param_%s_sd", p)]] <<- score_sd_wl[positions]
    data[[sprintf("Phact_wol_param_%s_sd", p)]] <<- score_sd_wol[positions]
    
    data[[sprintf("Phact_wl_param_%s_avg", p)]] <<- score_avg_wl[positions]
    data[[sprintf("Phact_wol_param_%s_avg", p)]] <<- score_avg_wol[positions]
  })
  
  ss <- mapply(function(ik){rep(scores[[sprintf("%s_num_of_trims", id)]][ik], 20)}, 1:length(scores[[sprintf("%s_num_of_trims", id)]]))
  data[["Num_of_Trims"]] <-  as.vector(ss) 
  ss <- mapply(function(ik){rep(scores[[sprintf("%s_num_of_diffaa", id)]][ik], 20)}, 1:length(scores[[sprintf("%s_num_of_trims", id)]]))
  data[["Num_of_Diff_AA"]] <-  as.vector(ss) 
  ss <- mapply(function(ik){rep(scores[[sprintf("%s_diffaa_trim", id)]][ik], 20)}, 1:length(scores[[sprintf("%s_num_of_trims", id)]]))
  data[["DiffAA_Trim"]] <- as.vector(ss) 
  
  ################MSA & Blosum Scores
  
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  data("BLOSUM62")
  blosum_mat <- BLOSUM62[amino_acids, amino_acids]
  blosum_mat <- rbind(blosum_mat, matrix(0, 1, 20))
  blosum_mat <- cbind(blosum_mat, matrix(0, 21, 1))
  colnames(blosum_mat)[21] <- "-"
  row.names(blosum_mat)[21] <- "-"
  data$Blosum62_score <- blosum_mat[cbind(data$Ref_AA, data$Alt_AA)]
  
  ############### ml features
  data$subs_acceptance <- sapply(1:length(positions), function(p){ml_features[[positions[p]]][["binary_accepted"]][[data$Alt_AA[[p]]]]})
  data$accepted_node_count <- sapply(1:length(positions), function(p){ml_features[[positions[p]]][["n_of_nodes_accepted"]][[data$Alt_AA[[p]]]]})
  data$dist_accepted <- sapply(1:length(positions), function(p){ml_features[[positions[p]]][["dist_accepted_min"]][[data$Alt_AA[[p]]]]})
  data$human_node_dist <- sapply(1:length(positions), function(p){ml_features[[positions[p]]][["human_node_dist"]]})
  
  human_node_probs_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["human_node_prob"]][data[p, "Ref_AA"]]}))
  colnames(human_node_probs_ref) <- "human_node_prob_ref"
  data <- cbind(data, human_node_probs_ref)
  
  human_node_probs_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["human_node_prob"]][data[p, "Alt_AA"]]}))
  colnames(human_node_probs_alt) <- "human_node_prob_alt"
  data <- cbind(data, human_node_probs_alt)
  
  data$delta_h_node_prob <- human_node_probs_ref - human_node_probs_alt
  
  el_i <- c()
  for (p in unique(positions)){
    if (msa[grep(id, row.names(msa)), p]=="-"){
      el_i <- c(el_i, which(data$Positions==p))
    }
  }
  positions <- positions[-el_i]
  data <- data[-el_i,]
  
  weighted_node_probs_ref <- t(as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["weighted_probs"]][,data[p, "Ref_AA"]]})))
  colnames(weighted_node_probs_ref) <- paste0("wp_", c(1:20)*5, "_ref")
  data <- cbind(data, weighted_node_probs_ref)
  
  weighted_node_probs_alt <- t(as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["weighted_probs"]][,data[p, "Alt_AA"]]})))
  colnames(weighted_node_probs_alt) <- paste0("wp_", c(1:20)*5, "_alt")
  data <- cbind(data, weighted_node_probs_alt)
  
  delta_w_node_probs <- weighted_node_probs_ref - weighted_node_probs_alt
  colnames(delta_w_node_probs) <- paste0("delta_wp_", c(1:20)*5)
  data <- cbind(data, delta_w_node_probs)
  
  
  tree_kernel <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["tree_kernel"]][data[p, "Alt_AA"]] }))
  colnames(tree_kernel) <- "tree_kernel"
  data <- cbind(data, tree_kernel)
  
  aa_class_ref <- t(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["aa_class"]]}))
  aa_class_alt <- aaComp(data$Alt_AA)
  aa_class_alt <- t(sapply(1:length(positions), function(i){aa_class_alt[[i]][,1]}))
  
  aa_class_comb <- aa_class_ref + aa_class_alt
  
  colnames(aa_class_ref) <- paste0(colnames(aa_class_ref), "_ref")
  colnames(aa_class_alt) <- paste0(colnames(aa_class_alt), "_alt")
  colnames(aa_class_comb) <- paste0(colnames(aa_class_comb), "_comb")
  
  data <- cbind(data, aa_class_ref, aa_class_alt, aa_class_comb)
  
  pos_freq_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["pos_freq"]][aa_to_num(data[p, "Ref_AA"])]}))
  pos_freq_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["pos_freq"]][aa_to_num(data[p, "Alt_AA"])]}))
  delta_pos_freq <- pos_freq_ref - pos_freq_alt
  
  #print(cbind(1:length(positions), pos_freq_ref))
  gap_freq <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["gap_freq"]]}))
  
  weighted_pos_freq_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["weighted_pos_freq"]][data[p, "Ref_AA"]]}))
  weighted_pos_freq_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["weighted_pos_freq"]][data[p, "Alt_AA"]]}))
  delta_weighted_pos_freq <- weighted_pos_freq_ref - weighted_pos_freq_alt
  
  
  data <- cbind(data, pos_freq_ref, pos_freq_alt, delta_pos_freq, gap_freq, weighted_pos_freq_ref, weighted_pos_freq_alt, delta_weighted_pos_freq)
  
  ind_increment_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["independent_increments_total"]][1,data[p, "Ref_AA"]]}))
  ind_increment_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["independent_increments_total"]][1,data[p, "Alt_AA"]]}))
  delta_ind_increments <- ind_increment_ref - ind_increment_alt
  
  data <- cbind(data, ind_increment_ref, ind_increment_alt, delta_ind_increments)
  
  accepted_ratio_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["accepted_ratio"]][[data[p, "Ref_AA"]]]}))
  accepted_ratio_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["accepted_ratio"]][[data[p, "Alt_AA"]]]}))
  data <- cbind(data, accepted_ratio_ref, accepted_ratio_alt)
  
  first_observed_node_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["first_observed_node"]][[data[p, "Ref_AA"]]]}))
  first_observed_node_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["first_observed_node"]][[data[p, "Alt_AA"]]]}))
  
  first_observed_node_alt[is.infinite(first_observed_node_alt)] <- 10000
  data <- cbind(data, first_observed_node_ref, first_observed_node_alt)
  
  
  allProbs_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["allProb_wl"]][1,data[p, "Ref_AA"]]}))
  allProbs_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["allProb_wl"]][1,data[p, "Alt_AA"]]}))
  delta_allProbs <- allProbs_ref - allProbs_alt
  
  data <- cbind(data, allProbs_ref, allProbs_alt, delta_allProbs)
  
  allDiff_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["allDiff_wl"]][1,data[p, "Ref_AA"]]}))
  allDiff_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["allDiff_wl"]][1,data[p, "Alt_AA"]]}))
  delta_allDiff <- allDiff_ref - allDiff_alt
  
  data <- cbind(data, allDiff_ref, allDiff_alt, delta_allDiff)
  
  score_aaprop_ref <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["score_aaprop"]][1,aa_to_num(data[p, "Ref_AA"])]}))
  score_aaprop_alt <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["score_aaprop"]][1,aa_to_num(data[p, "Alt_AA"])]}))
  delta_aaprop <- score_aaprop_ref - score_aaprop_alt
  
  score_aadiv <- as.matrix(sapply(1:length(positions), function(p){ml_features[[positions[p]]][["score_aadiv"]]}))
  
  data <- cbind(data, score_aaprop_ref, score_aaprop_alt, delta_aaprop, score_aadiv)
  
  ####################### Protein Level Features
  
  data$total_tree_length <- protein_level_features$total_tree_length
  data$optimal_loglikelihood <- protein_level_features$optimal_loglikelihood
  
  data$rate_of_no_ind_increment_ref <- protein_level_features$rate_of_no_ind_increment[data$Ref_AA]
  data$rate_of_no_ind_increment_alt <- protein_level_features$rate_of_no_ind_increment[data$Alt_AA]
  
  data$num_ancestors <- protein_level_features$num_ancestors
  data$mean_tree_length <- protein_level_features$mean_tree_length
  data$sd_tree_length <- protein_level_features$sd_tree_length
  
  data$total_gap_freq <- protein_level_features$total_gap_freq
  
  site_prop_rates <- matrix(protein_level_features$site_prop_rates, nrow = nrow(data), ncol = length(protein_level_features$site_prop_rates), byrow = T)
  colnames(site_prop_rates) <- names(protein_level_features$site_prop_rates)
  data <- cbind(data, site_prop_rates)
  
  data$overall_div <- as.numeric(protein_level_features[8])
  data$aver_evol_rate <- as.numeric(protein_level_features[10])
  data$total_passed <- as.numeric(protein_level_features[11])
  data$total_iter <- as.numeric(protein_level_features[12])
  data$num_of_const_site <- as.numeric(protein_level_features[13])
  data$time_wc <- as.numeric(protein_level_features[14])
  data$num_of_invariant_sites <- as.numeric(protein_level_features[15])
  data$num_of_pars_inf_sites <- as.numeric(protein_level_features[16])
  data$num_of_dist_sites <- as.numeric(protein_level_features[17])
  data$uncosnt_lh <- as.numeric(protein_level_features[18])
  data$num_free_params <- as.numeric(protein_level_features[19])
  data$AIC_score <- as.numeric(protein_level_features[20])
  data$AICc_score <- as.numeric(protein_level_features[21])
  data$BIC_score <- as.numeric(protein_level_features[22])
  data$internal_br_len <- as.numeric(protein_level_features[23])
  data$ratio_seq_pos <- as.numeric(protein_level_features[24])
  
  
  ######################### Protscale Scores
  
  scales <- names(protscale_scores$scores_aver)
  
  ref_aa <- match(data$Ref_AA, colnames(protscale_scores$scores_aver[[1]]))
  protscales_aver_ref <- list()
  protscales_tol_ref <- list()
  protscales_wo_phact_ref <- list()
  
  for(scale in scales) {
    protscales_aver_ref[[scale]] <- protscale_scores$scores_aver[[scale]][cbind(positions, ref_aa)]
    protscales_tol_ref[[scale]] <- protscale_scores$scores_tol[[scale]][cbind(positions, ref_aa)]
    protscales_wo_phact_ref[[scale]] <- protscale_scores$scores_wo_phact[[scale]][cbind(positions, ref_aa)]
  }
  protscales_aver_ref <- matrix(unlist(protscales_aver_ref), nrow = length(positions), ncol = length(protscales_aver_ref))
  protscales_tol_ref <- matrix(unlist(protscales_tol_ref), nrow = length(positions), ncol = length(protscales_tol_ref))
  protscales_wo_phact_ref <- matrix(unlist(protscales_wo_phact_ref), nrow = length(positions), ncol = length(protscales_wo_phact_ref))
  
  colnames(protscales_aver_ref) <- paste0(scales, "_aver_ref")
  colnames(protscales_tol_ref) <- paste0(scales, "_tol_ref")
  colnames(protscales_wo_phact_ref) <- paste0(scales, "_wo_phact_ref")
  
  data <- cbind(data, protscales_aver_ref, protscales_tol_ref, protscales_wo_phact_ref)
  
  alt_aa <- match(data$Alt_AA, colnames(protscale_scores$scores_aver[[1]]))
  protscales_aver_alt <- list()
  protscales_tol_alt <- list()
  protscales_wo_phact_alt <- list()
  
  for(scale in scales) {
    protscales_aver_alt[[scale]] <- protscale_scores$scores_aver[[scale]][cbind(positions, alt_aa)]
    protscales_tol_alt[[scale]] <- protscale_scores$scores_tol[[scale]][cbind(positions, alt_aa)]
    protscales_wo_phact_alt[[scale]] <- protscale_scores$scores_wo_phact[[scale]][cbind(positions, alt_aa)]
    
  }
  protscales_aver_alt <- matrix(unlist(protscales_aver_alt), nrow = length(positions), ncol = length(protscales_aver_alt))
  protscales_tol_alt <- matrix(unlist(protscales_tol_alt), nrow = length(positions), ncol = length(protscales_tol_alt))
  protscales_wo_phact_alt <- matrix(unlist(protscales_wo_phact_alt), nrow = length(positions), ncol = length(protscales_wo_phact_alt))
  
  colnames(protscales_aver_alt) <- paste0(scales, "_aver_alt")
  colnames(protscales_tol_alt) <- paste0(scales, "_tol_alt")
  colnames(protscales_wo_phact_alt) <- paste0(scales, "_wo_phact_alt")
  
  data <- cbind(data, protscales_aver_alt, protscales_tol_alt, protscales_wo_phact_alt)
  
  delta_protscales_aver <- protscales_aver_ref - protscales_aver_alt
  delta_protscales_tol <- protscales_tol_ref - protscales_tol_alt
  delta_protscales_wo_phact <- protscales_wo_phact_ref - protscales_wo_phact_alt
  
  colnames(delta_protscales_aver) <- paste0("delta_", scales, "_aver")
  colnames(delta_protscales_tol) <- paste0("delta_", scales, "_tol")
  colnames(delta_protscales_wo_phact) <- paste0("delta_", scales, "_wo_phact")
  
  data <- cbind(data, delta_protscales_aver, delta_protscales_tol, delta_protscales_wo_phact)
  
  ##############
  ##############
  
  ###############
  c_names <- colnames(data)
  data[is.infinite(data$first_observed_node_ref), "first_observed_node_ref"] <- 10000
  data[is.na(data$weighted_pos_freq_ref), "weighted_pos_freq_ref"] <- 0 
  data[is.na(data$weighted_pos_freq_alt), "weighted_pos_freq_alt"] <- 0 
  data[is.na(data$delta_weighted_pos_freq), "delta_weighted_pos_freq"] <- 0 
  
  for(param in parameters){
    data[[sprintf("Phact_wl_param_%s_delta", param)]] <- data[[sprintf("Phact_wl_param_%s_ref", param)]] - data[[sprintf("Phact_wl_param_%s_alt", param)]]
    data[[sprintf("Phact_wol_param_%s_delta", param)]] <- data[[sprintf("Phact_wol_param_%s_ref", param)]] - data[[sprintf("Phact_wol_param_%s_alt", param)]]
    
    data[[sprintf("Phact_wl_param_%s_ratio", param)]] <- data[[sprintf("Phact_wl_param_%s_alt", param)]] / data[[sprintf("Phact_wl_param_%s_ref", param)]]
    data[[sprintf("Phact_wol_param_%s_ratio", param)]] <- data[[sprintf("Phact_wol_param_%s_alt", param)]] / data[[sprintf("Phact_wol_param_%s_ref", param)]]
  }
  
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  load("./aa_scales.RData")
  aa_scales <- as.data.frame(aa_scales)
  rownames(aa_scales) <- aa_scales[,1]
  aa_scales_ref <- (aa_scales[data$Ref_AA, -1])
  aa_scales_alt <- (aa_scales[data$Alt_AA, -1])
  aa_scales_delta <- aa_scales_ref - aa_scales_alt
  colnames(aa_scales_ref) <- paste0(colnames(aa_scales_ref), "_ref")
  colnames(aa_scales_alt) <- paste0(colnames(aa_scales_alt), "_alt")
  colnames(aa_scales_delta) <- paste0(colnames(aa_scales_delta), "_delta")
  
  data <- cbind(data, aa_scales_ref, aa_scales_alt, aa_scales_delta)
  #data$first_position <- 1*(data$Ref_AA == "M" & data$Position == 1)
  
  features <- colnames(data)
  features <- gsub(",", "_", features)
  colnames(data) <- features
  
  data[,"accepted_node_count"] <- log(data[,"accepted_node_count"] + 1e-100)
  
  ##############
  ##############
  if (length(keep_gaps)>=1){
    for (i in 1:length(keep_gaps)){
      mat <- data[1:20,]
      mat[,2] <- keep_gaps[i]
      mat[,3] <- "-"
      mat[,5:length(mat[1,])] <- 0
      
      if (keep_gaps[i]==1){
        data[,2] <- as.numeric(data[,2]) + 1
        data <- rbind(mat, data)
      } else {
        k <- (keep_gaps[i]-1)*20
        part1 <- data[1:k,]
        part2 <- data[(k+1):length(data[,1]),]
        part2[,2] <- as.numeric(part2[,2]) + 1
        data <- rbind(part1, mat, part2)
      }
    }
  }
  
  
  save("data", file = sprintf("%s/%s.RData", save_path, id))
}




