#!/usr/bin/env Rscript
library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(Peptides)
source("./compute_PHACT.R")
source("./compute_input_features.R")
source("./compute_weight.R")

args = commandArgs(trailingOnly=TRUE)

uniprot_id <- args[6]
save_path <- args[9]

if(dir.exists(sprintf("%s", save_path)) == FALSE) {
  dir.create(sprintf("%s",save_path))
}

if(dir.exists(sprintf("%s/%s", save_path, uniprot_id)) == FALSE) {
  dir.create(sprintf("%s/%s",save_path, uniprot_id))
}

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

compute_score <- function(file_nwk, file_rst, file_fasta, file_log, file_iqtree, output_name, human_id, pos_chosen, parameters, aa_scales_file) {
  
  # Read tree file
  tr_org <- read.tree(file_nwk)
  x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
  x[,1] <- str_remove(x[,1], "Node")
  colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )
  # Tree_info: node-node, node-leaf connections
  tree_info <- as.data.frame(as_tibble(tr_org))
  
  # Read fasta file, MSA
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  
  # connections_1: Parent node, connections_2: connected node/leaf
  connections_1 <- tree_info$parent
  connections_2 <- tree_info$node
  
  # Names of leaves
  names_all <- tr_org[["tip.label"]]
  msa <- msa[names_all, ]
  # Number of total leaves&nodes
  num_leaves <- length(tr_org[["tip.label"]])
  num_nodes <- tr_org[["Nnode"]]
  
  # Distance between leaves & nodes
  dd_node <- dist.nodes(tr_org)
  dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
  dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
  
  # Human position (leaf & node)
  h_name <- human_id
  human_codeml <- names_all[grep(pattern = h_name, x = names_all, fixed = TRUE)]
  leaf_human <- tree_info[which(tree_info$label == human_codeml), "node"]
  human_plc <- leaf_human
  node_human <- tree_info[which(tree_info$label == human_codeml), "parent"]
  nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
  names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
  
  # Total number of positions from ancestralProbs file
  total_pos <- max(x$Site)
  
  # Chosen positions (all or some)
  if (pos_chosen[1] == "all"){
    positions <- 1:total_pos
    score_all <- matrix(0, total_pos, 21)
  } else {
    positions <- pos_chosen
    score_all <- matrix(0, length(positions), 21)
  }
  
  human_gaps <- which(msa[match(human_codeml, row.names(msa))]=="-")
  positions <- setdiff(positions, human_gaps)
  
    
  ####################################################
  ####################################################
  
  # Connections between leaves & nodes
  chosen_leaves <- tree_info[1:num_leaves,c("parent", "node")]
  # Connections between nodes & nodes
  chosen_nodes <- tree_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
  leaf_names <- tree_info$label
  
  human_leaf_len <- as.double(tree_info[human_plc, "branch.length"])
  
  if (num_nodes == 1) {
    d_n <- dist_node + human_leaf_len
  } else {
    d_n <- dist_node[as.character(node_human),] + human_leaf_len
  }
  
  d_l <- dist_leaf[leaf_human,]
  
  # chosen_nodes2: ordered connections (for probability differences)
  chosen_nodes2 <- matrix(0, num_nodes-1, 2)
  
  n1 <- as.numeric(chosen_nodes$parent)
  n2 <- as.numeric(chosen_nodes$node)
  dist_f <- d_n[as.character(n1)]
  dist_s <- d_n[as.character(n2)]
  
  # chosen_nodes2: ordered connections (for probability differences)
  chosen_nodes2[which(dist_f < dist_s), 1] <- n2[which(dist_f < dist_s)]
  chosen_nodes2[which(dist_f < dist_s), 2] <- n1[which(dist_f < dist_s)]
  
  chosen_nodes2[which(dist_f >= dist_s), 1] <- n1[which(dist_f >= dist_s)]
  chosen_nodes2[which(dist_f >= dist_s), 2] <- n2[which(dist_f >= dist_s)]
  
  # Number of nodes between nodes & leaf of human
  nodes_conn <- numeric(num_nodes)
  nodes_conn[node_human-num_leaves] <- 1
  names(nodes_conn) <- names(d_n)
  chs <- c()
  chs2 <- c()
  ##########################
  inds <- chosen_nodes2[chosen_nodes2[,2]==(node_human),1]
  nodes_conn[as.character(inds)] <- 2
  chs <- inds
  
  s0 <- sapply(3:num_leaves, function(i){
    for (j in chs){
      inds <- chosen_nodes2[chosen_nodes2[,2]==j,1]
      if (length(inds)!=0){
        nodes_conn[as.character(inds)] <<- i
        chs2 <- c(chs2, inds)
      }
    }
    chs <<- chs2
    chs2 <- c()
  })
  
  # Number of nodes between leaves & leaf of human
  if (num_nodes == 1) {
    leaves_conn <- nodes_conn*matrix(1, 1, num_leaves)
  } else {
    leaves_conn <- nodes_conn[as.character(chosen_leaves[,1])]
  }  
  #########################################
  
  #########aa_scales
  load(aa_scales_file)
  #########
  
  polarity_polar <- c("C", "D", "E", "H", "K", "N", "Q", "R", "S", "T", "Y")
  polarity_nonpolar <- c("A", "F", "G", "I", "L", "M", "P", "V", "W")
  
  charge_positive <- c("H", "K", "R")
  charge_negative <- c("D", "E")
  charge_neutral <- c("A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", "V", "W", "Y")
  
  aromaticity_aliphatic <- c("I", "L", "V")
  aromaticity_aromatic <- c("F", "H", "W", "Y")
  aromaticity_neutral <- c("A", "C", "D", "E", "G", "K", "M", "N", "P", "Q", "R", "S", "T")
  
  size_small <- c("A", "G", "P", "S")
  size_medium <- c("D", "N", "T")
  size_large <- c("C", "E", "F", "H", "I", "K", "L", "M", "Q", "R", "V", "W", "Y")
  
  electronic_strongdonor <- c("A", "D", "E", "P")
  electronic_weak_donor <- c("I", "L", "V")
  electronic_neutral <- c("C", "G", "H", "S", "W")
  electronic_weak_acceptor <- c("F", "M", "Q", "T", "Y")
  electronic_strong_acceptor <- c("K", "N", "R")
  
  type_aas <- matrix(0, 5, 5)
  type_aas[1,1] <- toString(polarity_polar)
  type_aas[1,2] <- toString(polarity_nonpolar)
  
  type_aas[2,1] <- toString(charge_positive)
  type_aas[2,2] <- toString(charge_negative)
  type_aas[2,3] <- toString(charge_neutral)
  
  type_aas[3,1] <- toString(aromaticity_aliphatic)
  type_aas[3,2] <- toString(aromaticity_aromatic)
  type_aas[3,3] <- toString(aromaticity_neutral)
  
  type_aas[4,1] <- toString(size_small)
  type_aas[4,2] <- toString(size_medium)
  type_aas[4,3] <- toString(size_large)
  
  type_aas[5,1] <- toString(electronic_strongdonor)
  type_aas[5,2] <- toString(electronic_weak_donor)
  type_aas[5,3] <- toString(electronic_neutral)
  type_aas[5,4] <- toString(electronic_weak_acceptor)
  type_aas[5,5] <- toString(electronic_strong_acceptor)
  
  ##########
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  
  score_norm <- t(mapply(function(ps){position_score(ps, x, msa, human_id, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes, nodes_raxml, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, parameters)}, rep(positions)))
  
  pl <- 1
  scores <- list()
  for (p in 1:length(parameters)) {
    score_norm_with_leaf <- matrix(unlist(score_norm[ ,pl]), nrow = length(positions), ncol = 20, byrow = TRUE)
    score_norm_without_leaf <- matrix(unlist(score_norm[ ,(pl+1)]), nrow = length(positions), ncol = 20, byrow = TRUE)
    score_norm_diversity <- matrix(unlist(score_norm[ ,(pl+2)]), nrow = length(positions), ncol = 1, byrow = TRUE)
    score_norm_nogap <- matrix(unlist(score_norm[ ,(pl+3)]), nrow = length(positions), ncol = 20, byrow = TRUE)
    
    pl <- pl + 4
    
    score_norm_with_leaf <- cbind(1:length(positions), score_norm_with_leaf)
    score_norm_without_leaf <- cbind(1:length(positions), score_norm_without_leaf)
    score_norm_nogap <- cbind(1:length(positions), score_norm_nogap)
    
    colnames(score_norm_with_leaf) <- c("Pos/AA", num_to_aa(1:20))
    colnames(score_norm_without_leaf) <- c("Pos/AA", num_to_aa(1:20))
    colnames(score_norm_nogap) <- c("Pos/AA", num_to_aa(1:20))
    
    filename <- ifelse(parameters[p] == "0", "max05", ifelse(parameters[p] == "X", "max05_Gauss", parameters[p]))
    filename <- ifelse(parameters[p] == "X", "max05_Gauss", parameters[p])
    scores[[sprintf("%s_wl_param_%s", output_name, filename)]] <- score_norm_with_leaf
    scores[[sprintf("%s_wol_param_%s", output_name, filename)]] <- score_norm_without_leaf
    scores[[sprintf("%s_diversity_%s", output_name, filename)]] <- score_norm_diversity
    scores[[sprintf("%s_nogap_param_%s", output_name, filename)]] <- score_norm_nogap
  }
  trims <- matrix(unlist(score_norm[ ,(pl)]), nrow = length(positions), ncol = 1, byrow = TRUE)
  pl <- pl + 1
  diversity_aa <- matrix(unlist(score_norm[ ,(pl)]), nrow = length(positions), ncol = 1, byrow = TRUE)
  pl <- pl + 1
  diversity_aa2 <- matrix(unlist(score_norm[ ,(pl)]), nrow = length(positions), ncol = 1, byrow = TRUE)
  
  scores[[sprintf("%s_num_of_trims", output_name)]] <- trims
  scores[[sprintf("%s_num_of_diffaa", output_name)]] <- diversity_aa
  scores[[sprintf("%s_diffaa_trim", output_name)]] <- diversity_aa2
  print("saving scores") 
  save("scores", file = sprintf("%s/%s/%s_scores.RData", save_path, output_name, output_name))
  
  ml_features_prot_scale <- sapply(positions, function(ps){calculate_ml_features(ps, x, msa, type_aas, human_codeml, tr_org, tree_info, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, chosen_leaves, chosen_nodes2, d_n, d_l, aa_scales, nodes_conn, leaves_conn)})
  ml_features <- ml_features_prot_scale[1,]
  protscale_scores <- ml_features_prot_scale[2,]
  scale_scores_wo_phact <- ml_features_prot_scale[3,]
  protscale_scores_aver <- list()
  protscale_scores_tol <- list()
  protscale_scores_wo_phact <- list()
  for (i in names(protscale_scores[[1]])){
    protscale_scores_aver[[i]] <- t(sapply(1:length(protscale_scores), function(j){protscale_scores[[j]][[i]]$scale_average}))
    protscale_scores_tol[[i]] <- t(sapply(1:length(protscale_scores), function(j){protscale_scores[[j]][[i]]$scale_tol}))
    protscale_scores_wo_phact[[i]] <- t(sapply(1:length(protscale_scores), function(j){scale_scores_wo_phact[[j]][[i]]}))
    
    colnames(protscale_scores_aver[[i]]) <- num_to_aa(1:20)
    colnames(protscale_scores_tol[[i]]) <- num_to_aa(1:20)
    colnames(protscale_scores_wo_phact[[i]]) <- num_to_aa(1:20)
    
  }
  protscale_scores <- list(scores_aver = protscale_scores_aver, scores_tol = protscale_scores_tol, scores_wo_phact = protscale_scores_wo_phact)
  
  log_file <- readLines(file_log)
  
  total_tree_length <- log_file[grep("Total tree length", log_file)]
  total_tree_length <- as.numeric(str_split(total_tree_length, ":", simplify = T)[2])
  
  optimal_loglikelihood <- log_file[grep("Optimal log-likelihood", log_file)]
  optimal_loglikelihood <- as.numeric(str_split(optimal_loglikelihood, ":", simplify = T)[2])
  
  site_prop_rates <- log_file[grep("Site proportion and rates", log_file)]
  site_prop_rates <- regmatches(site_prop_rates, gregexpr("(?<=\\().*?(?=\\))", site_prop_rates, perl=T))[[1]]
  site_prop_rates <- as.numeric(sapply(site_prop_rates, function(x){str_split(x,",",simplify = T)}))
  names(site_prop_rates) <- c("prop1", "rate1", "prop2", "rate2", "prop3", "rate3", "prop4", "rate4")
  
  aver_evol_rate <- sum(site_prop_rates[c(1,3,5,7)]*site_prop_rates[c(2,4,6,8)])
  
  total_passed <- log_file[grep("TOTAL", log_file)][1]
  total_passed <- str_split(total_passed, "L", simplify = T)[2]
  total_passed <- as.numeric(substr(total_passed, 1, (which(strsplit(total_passed, "")[[1]]=="%")-1)))
  
  total_iter <- log_file[grep("Total number of iterations", log_file)]
  total_iter <- as.numeric(str_split(total_iter, ":", simplify = T)[2])
  
  iqtree_file <- readLines(file_iqtree)
  num_of_const_sites <- iqtree_file[grep("Number of constant sites", iqtree_file)]
  num_of_const_sites <- str_split(num_of_const_sites, ":", simplify = T)[2]
  num_of_const_sites <- as.numeric(substr(num_of_const_sites, 1, (which(strsplit(num_of_const_sites, "")[[1]]=="(")-1)))
  
  time_wc <- iqtree_file[grep("Total wall-clock time used", iqtree_file)]
  time_wc <- str_split(time_wc, ":", simplify = T)[2]
  time_wc <- as.numeric(substr(time_wc, 1, (which(strsplit(time_wc, "")[[1]]=="s")-1)))
  
  num_of_invariant_sites <- iqtree_file[grep("Number of invariant", iqtree_file)]
  num_of_invariant_sites <- str_split(num_of_invariant_sites, ":", simplify = T)[2]
  num_of_invariant_sites <- as.numeric(substr(num_of_invariant_sites, 1, (which(strsplit(num_of_invariant_sites, "")[[1]]=="(")-1)))
  
  num_of_pars_inf_sites <- iqtree_file[grep("Number of parsimony informative sites", iqtree_file)]
  num_of_pars_inf_sites <- as.numeric(str_split(num_of_pars_inf_sites, ":", simplify = T)[2])
  
  num_of_dist_sites <- iqtree_file[grep("Number of distinct site patterns", iqtree_file)]
  num_of_dist_sites <- as.numeric(str_split(num_of_dist_sites, ":", simplify = T)[2])
  
  uncosnt_lh <- iqtree_file[grep("Unconstrained log-likelihood", iqtree_file)]
  uncosnt_lh <- as.numeric(str_split(uncosnt_lh, ":", simplify = T)[2])
  
  num_free_params <- iqtree_file[grep("Number of free parameters", iqtree_file)]
  num_free_params <- as.numeric(str_split(num_free_params, ":", simplify = T)[2])
  
  AIC_score <- iqtree_file[grep("Akaike information criterion", iqtree_file)][1]
  AIC_score <- as.numeric(str_split(AIC_score, ":", simplify = T)[2])
  AICc_score <- iqtree_file[grep("Corrected Akaike information criterion", iqtree_file)]
  AICc_score <- as.numeric(str_split(AICc_score, ":", simplify = T)[2])
  BIC_score <- iqtree_file[grep("Bayesian information criterion", iqtree_file)]
  BIC_score <- as.numeric(str_split(BIC_score, ":", simplify = T)[2])
  
  internal_br_len <- iqtree_file[grep("Sum of internal branch lengths", iqtree_file)]
  internal_br_len <- str_split(internal_br_len, ":", simplify = T)[2]
  internal_br_len <- as.numeric(substr(internal_br_len, 1, (which(strsplit(internal_br_len, "")[[1]]=="(")-1)))
  
  ratio_seq_pos <- iqtree_file[grep("Input data: ", iqtree_file)]
  ratio_seq_pos <- substr(ratio_seq_pos, 12, (nchar(ratio_seq_pos)-17))
  ratio_seq_pos <- as.numeric(substr(ratio_seq_pos, 1, (which(strsplit(ratio_seq_pos, "")[[1]]=="s")-1)))/
    as.numeric(substr(ratio_seq_pos, (which(strsplit(ratio_seq_pos, "")[[1]]=="h")+1), nchar(ratio_seq_pos)))
  
  rate_of_no_ind_increment <- colSums(t(sapply(1:length(positions), function(p){1*(ml_features[[p]]$independent_increments_total == 0)}))) / max(positions)
  names(rate_of_no_ind_increment) <- num_to_aa(1:20)
  
  mean_tree_length <- mean(tree_info$branch.length, na.rm = T)
  sd_tree_length <- sd(tree_info$branch.length, na.rm = T)
  
  total_gap_freq <- sum(msa == "-") / (max(positions)*num_nodes)
  
  ind_incs_all <- matrix(0, length(positions), 20)
  for (ps in 1:length(positions)){
    ind_incs_all[ps, ] <- ml_features[[ps]]$independent_increments_total
  }
  ind_inc_tol <- sum(ind_incs_all)/length(positions)
  
  protein_level_features <- list(total_tree_length = total_tree_length, optimal_loglikelihood = optimal_loglikelihood, 
                                 rate_of_no_ind_increment = rate_of_no_ind_increment, num_ancestors = num_nodes, mean_tree_length = mean_tree_length,
                                 sd_tree_length = sd_tree_length, total_gap_freq = total_gap_freq, overall_div = ind_inc_tol,
                                 site_prop_rates = site_prop_rates, aver_evol_rate = aver_evol_rate,
                                 total_passed = total_passed, total_iter = total_iter, num_of_const_sites = num_of_const_sites,
                                 time_wc = time_wc, num_of_invariant_sites = num_of_invariant_sites,
                                 num_of_pars_inf_sites = num_of_pars_inf_sites, num_of_dist_sites = num_of_dist_sites,
                                 uncosnt_lh = uncosnt_lh, num_free_params = num_free_params, 
                                 AIC_score = AIC_score, AICc_score = AICc_score, BIC_score = BIC_score, 
                                 internal_br_len = internal_br_len, ratio_seq_pos = ratio_seq_pos)
  print("saving files")
  save("ml_features", file = sprintf("%s/%s/%s_ml_features.RData", save_path, output_name, output_name))
  save("protscale_scores", file = sprintf("%s/%s/%s_protscale_scores.RData", save_path, output_name, output_name))

  save("protein_level_features", file = sprintf("%s/%s/%s_protein_level_features.RData", save_path, output_name, output_name))
}


compute_score(file_nwk=args[1],file_rst=args[2],file_fasta=args[3],file_log=args[4],file_iqtree=args[5], output_name=args[6],human_id=args[7],'all', parameters = args[8], aa_scales_file=args[10])

