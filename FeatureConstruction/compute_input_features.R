calculate_ml_features <- function(ps, version, x, msa, type_aas, human_codeml, tr_org, tree_info, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, chosen_leaves, chosen_nodes2, d_n, d_l, aa_scales, nodes_conn, leaves_conn) {
  position <- ps
  
  version <- 0
  #version <- as.numeric(version)
  
  aa_scales_org <- aa_scales
  if (version==1){
    aa_scales <- aa_scales[match(num_to_aa(1:20), unlist(aa_scales[,1])),]
  } else {
    aa_scales <- aa_scales[match(num_to_aa(1:20), unlist(aa_scales[,1])),]
    aa_scales[,2:27] <- scale(aa_scales[,2:27])
  }
  
  num_nodes_prev <- num_nodes
  num_leaves_prev <- num_leaves
  nodes_raxml_prev <- nodes_raxml
  
  b1 <- position + total_pos*(0:(num_nodes-1))
  TT <- x[b1,]
  
  node_info <- as.numeric(TT[,1])
  sort_node_info <- sort(node_info, decreasing = F, index.return=T)
  TT <- TT[sort_node_info$ix,]
  
  matrix_prob <- matrix(0, num_nodes, 20)
  
  probs <- data.matrix((TT[, (4:ncol(TT))]))
  rownames(probs) <- NULL
  rr <- aa_to_num(colnames(x)[4:ncol(TT)])
  matrix_prob[,rr] <- probs
  matrix_prob <- matrix_prob[nodes_raxml,]
  if (length(matrix_prob)==20){
    names(matrix_prob) <- names(sort(rr))
  } else {
    colnames(matrix_prob) <- names(sort(rr))
  }
  
  msa_upd <- msa
  names_msa <- rownames(msa)
  trims <- names_msa[which(msa[,ps]=="-")]
  trims <- unique(trims)

  if (is.element(human_codeml, trims)){
    trims <- trims[-which(trims==human_codeml)]
  }
  
  if (length(trims)>=1) {
    tree_new <- drop.tip(tr_org, trims)
    tree_new_info <- as.data.frame(as_tibble(tree_new))
    
    num_leaves <- length(tree_new[["tip.label"]])
    num_nodes <- tree_new[["Nnode"]]
    
    nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_new_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
    names(nodes_raxml) <- tree_new_info[num_leaves+1:num_nodes, "node"]
    
    elim_nodes <- setdiff(nodes_raxml_prev, nodes_raxml)
    if (length(elim_nodes)>0){
      els <- t(mapply(function(jj){which(as.numeric(nodes_raxml_prev) == elim_nodes[jj])}, rep(1:length(elim_nodes))))
      matrix_prob <- matrix_prob[-els,]
    }
    
    if (length(trims)>=1){
      msa_upd <- msa_upd[-t(mapply(function(jj){which(rownames(msa_upd)==trims[jj])}, rep(1:length(trims)))),]
    }
    
    # Human position (leaf & node)
    leaf_human <- tree_new_info[which(tree_new_info$label == human_codeml), "node"]
    human_plc <- leaf_human
    node_human <- tree_new_info[which(tree_new_info$label == human_codeml), "parent"]
    
    # Distance between leaves & nodes
    dd_node <- dist.nodes(tree_new)
    dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
    dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
    
    # Connections between leaves & nodes
    chosen_leaves <- tree_new_info[1:num_leaves,c("parent", "node")]
    chosen_nodes <- tree_new_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
    leaf_names <- tree_new_info$label
    
    human_leaf_len <- as.double(tree_new_info[human_plc, "branch.length"])
    
    if (num_nodes == 1) {
      if (num_leaves == 1) {
        d_n <- dist_node + human_leaf_len
        d_l <- dist_leaf
      } else {
        d_n <- dist_node + human_leaf_len
        d_l <- dist_leaf[leaf_human,]
      }
    } else {
      d_n <- dist_node[as.character(node_human),] + human_leaf_len
      d_l <- dist_leaf[leaf_human,]
    }
    
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
    
  } else {
    tree_new <- tr_org
    tree_new_info <-  tree_info
  }
  
  if (num_leaves>1){
    position_vec <- msa_upd[, ps]
    position_num <- aa_to_num(position_vec)
    aa_f <- position_num[human_plc]
    ref_aa <- msa_upd[human_plc, ps]
  } else {
    position_vec <- msa_upd[ps]
    position_num <- aa_to_num(position_vec)
    aa_f <- position_num[human_plc]
    ref_aa <- msa_upd[ps]
  }
  gaps <- which(position_num == 21)
  
  parameter <- "CountNodes_3"
  weights <- weight_fnc(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn, mxx) # Updated with related parameters
  weight_leaf <- weights[1:num_leaves]
  weight_node <- tail(weights,num_nodes)
  
  if (num_nodes == 1) {
    scales <- colnames(aa_scales)[2:(dim(aa_scales)[2])]
    scale_list <- list()
    scale_matrices <- list()
    scale_leaf_probs <- list()
    
    prob_leaves <- matrix(0, num_leaves, 20)
    prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
    
    for(scale in scales){
      scale_list[[scale]] <- aa_scales[[scale]][match(names(matrix_prob), aa_scales$aa_names)]
      scale_matrices[[scale]] <- matrix_prob * matrix(scale_list[[scale]], nrow = 1, ncol = 20, byrow = T)
      scale_leaf_probs[[scale]] <- prob_leaves * matrix(scale_list[[scale]], nrow = num_leaves, ncol = 20, byrow = T)
    }
    scale_scores_wo_phact <- lapply(scale_matrices, function(scale_mat){ 
      scores <- colSums(scale_mat * matrix(weight_node, nrow = nrow(scale_mat), ncol = ncol(scale_mat), byrow = F))
    })
  } else {
    scales <- colnames(aa_scales)[2:(dim(aa_scales)[2])]
    scale_list <- list()
    scale_matrices <- list()
    scale_leaf_probs <- list()
    
    prob_leaves <- matrix(0, num_leaves, 20)
    prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
    
    for(scale in scales){
      scale_list[[scale]] <- aa_scales[[scale]][match(colnames(matrix_prob), aa_scales$aa_names)]
      scale_matrices[[scale]] <- matrix_prob * matrix(scale_list[[scale]], nrow = nrow(matrix_prob), ncol = 20, byrow = T)
      scale_leaf_probs[[scale]] <- prob_leaves * matrix(scale_list[[scale]], nrow = num_leaves, ncol = 20, byrow = T)
    }
    scale_scores_wo_phact <- lapply(scale_matrices, function(scale_mat){ 
      scores <- colSums(scale_mat * matrix(weight_node, nrow = nrow(scale_mat), ncol = ncol(scale_mat), byrow = F)) #1/dn ve weight dene!!!
    })
  }
  
  if (num_nodes == 1) {
    diff_leaves <- matrix(0, num_leaves, 20)
    diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
    diff_nodes <- matrix(0, 1, 20)
    vect_human <- matrix_prob
    diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
  } else {
    diff_leaves <- matrix(0, num_leaves, 20)
    diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
    diff_nodes <- matrix(0, num_nodes-1, 20)
    diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
    vect_human <- matrix_prob[node_human - num_leaves,]
    diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
  }
  
  ################################
  ###############  INDEPENDENT INC 
  ind_inc_nodes <- matrix(0,1,20)
  colnames(ind_inc_nodes) <- num_to_aa(1:20)
  
  ind_inc_leaves <- matrix(0,1,20)
  colnames(ind_inc_leaves) <- num_to_aa(1:20)
  
  if (num_nodes>1) {
    if (length(diff_nodes)>20) {
      for (j in 1:20){
        plc <- c()
        for (i in 1:(num_nodes-1)){
          if (diff_nodes[i,j]>0.01) {
            plc <- c(plc, i)
          }
        }
        if (length(plc)>0){
          for (i in 1:length(plc)) {
            val <- diff_nodes[plc[i],j]
            lc <- plc[i]
            left <- chosen_nodes2[lc,1]
            right <- chosen_nodes2[lc,2]
            
            p2 <- which(chosen_nodes2[,1]==right)
            if (length(p2) !=0){
              check_val <- diff_nodes[p2,j]
              if ((max(check_val)<0 && val>0.01) || (max(check_val)<0.005 && val>0.99)){
                ind_inc_nodes[j] <- ind_inc_nodes[j]+1
              }
            }
          }
        }
      }
    }
  }
  
  if (num_leaves>1) {
    for (j in 1:20){
      plc <- c()
      for (i in 1:(num_leaves)){
        if (diff_leaves[i,j]>0.01) {
          plc <- c(plc, i)
        }
      }
      if (length(plc)>0){
        for (i in 1:length(plc)) {
          val <- diff_leaves[plc[i],j]
          if (chosen_leaves[plc[i],2]==human_plc){
            ind_inc_leaves[j] = ind_inc_leaves[j] + 1
            break
          }
          lc <- plc[i]
          left <- chosen_leaves[lc,2]
          right <- chosen_leaves[lc,1]
          
          p2 <- which(chosen_nodes2[,1]==right)
          if (length(p2)!=0){
            if (length(diff_nodes)>20){
              check_val <- diff_nodes[p2,j]
            } else {
              check_val <- diff_nodes[j]
            }
            if  ((max(check_val)<0 && val>0.01) || (max(check_val)<0.005 && val>0.99)) {
              ind_inc_leaves[j] <- ind_inc_leaves[j]+1
            }
          }
        }
      }
    }
  }
  
  if (aa_f!=21){
    ind_inc_leaves[aa_f] <- ind_inc_leaves[aa_f] + 1
  }
  ind_inc_total <- ind_inc_leaves + ind_inc_nodes
  colnames(ind_inc_total) <- num_to_aa(1:20)
  
  ################################
  ######## SCALE SCORES
  aa_scales <- aa_scales_org
  scale_scores <- lapply(scales, function(scale){
    scale_val <- aa_scales[[scale]]
    # mat_prob <- scale_matrices[[scale]]
    # mat_prob <- matrix(unlist(mat_prob), num_nodes, 20)
    # prob_leaves_scaled <- scale_leaf_probs[[scale]]
    diff_leaves <- matrix(0, num_leaves, 20)
    if (num_nodes == 1) {
      diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
      diff_nodes <- matrix(0, 1, 20)
      vect_human <- matrix_prob
    } else {
      diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
      diff_nodes <- matrix(0, num_nodes-1, 20)
      diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
      vect_human <- matrix_prob[node_human - num_leaves,]
    }
    diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
    
    ####################
    aa_f_val <- 1
    
    score <- matrix(0,1,20)
    if (num_leaves == 1) {
      if (aa_f<=20){
        score[aa_f] <- aa_f_val
      }
      score_without_leaf <- score
    } else {
      s1 <- sapply(1:20, function(ii){
        if (length(diff_nodes)==20){
          dif_pr <- matrix(0,1,20)
        } else {
          dif_pr <- diff_nodes[1:(num_nodes-1),ii]
          dif_pr[dif_pr<0] <- 0
          
          sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
          score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr)
        }
      })
      
      if (length(diff_nodes)==20 && num_nodes==1){
        vect_human <- matrix_prob
      } else {
        vect_human <- matrix_prob[node_human - num_leaves,]
      }
      
      if (aa_f != 21) {
        vect_human[aa_f]<-0
      }
      score_without_leaf <- score
      score_without_leaf <- score_without_leaf + weight_node[(node_human-num_leaves)]*vect_human
      score <- score + weight_node[(node_human-num_leaves)]*vect_human
      
      s2 <- sapply(1:20, function(ii){
        diff_lf <- diff_leaves[1:num_leaves,ii]
        diff_lf[diff_lf<0] <- 0
        
        s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
        score[ii] <<- score[ii] + s1
      })
      
      if (aa_f != 21){
        score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*aa_f_val
        score_without_leaf[aa_f] <- score_without_leaf[aa_f] + weight_leaf[human_plc]*aa_f_val
      }
    }
    
    score_norm <- sum(score/(num_leaves+num_nodes)*scale_val)/sum(score/(num_leaves+num_nodes))*matrix(1,1,20)
    
    scale_tolerance <- 0
    prob_leaves <- prob_leaves*t(matrix(rep(as.numeric(scale_val), num_leaves), 20, num_leaves))
    matrix_prob <- matrix_prob*t(matrix(rep(as.numeric(scale_val), num_nodes), 20, num_nodes))
    diff_leaves <- matrix(0, num_leaves, 20)
    if (num_nodes == 1) {
      diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
      diff_nodes <- matrix(0, 1, 20)
      vect_human <- matrix_prob
    } else {
      diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
      diff_nodes <- matrix(0, num_nodes-1, 20)
      diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
      vect_human <- matrix_prob[node_human - num_leaves,]
    }
    diff_leaves[human_plc,] <- matrix(0,1,20)
    
    if (length(diff_nodes)==20){
      diff_nodes <- sum(diff_nodes)
    } else {
      diff_nodes <- rowSums(diff_nodes)
    }
    diff_leaves <- rowSums(diff_leaves)
    
    scale_tolerance <- 0
    
    if (num_nodes==1){
      sel_node <- 1
    } else {
      sel_node <- chosen_nodes2[1:length(diff_nodes), 1] - num_leaves
    }
    
    scale_tolerance_node <- sum(abs(weight_node[sel_node] * diff_nodes))
    scale_tolerance_leaf <- sum(abs(weight_leaf* diff_leaves))
    scale_tolerance <- (scale_tolerance_leaf + scale_tolerance_node)/(num_leaves+num_nodes)
    
    score_wol_norm <- scale_tolerance*matrix(1,1,20)
    scores <- list(scale_average = score_wol_norm, scale_tol = score_norm)
    
  })
  names(scale_scores) <- scales
  
  
  ########## POLARITY SCALE ###############################################
  
  score_aa_prop <- list()
  if (num_leaves == 1){
    score_aa <- matrix(0,1,20)
    aa_f <- aa_to_num(msa_upd[ps])
    if (aa_f!=21){
      score_aa[aa_f] <- 1
    }
    diversity <- 0
    position_vec <- msa_upd[ps]
    position_num <- aa_to_num(position_vec)
    
    pl <- 1
    score_aa_prop$score_aa <- score_aa/(num_nodes+num_leaves)
    score_aa_prop$div <- 0

  } else {
    
    position_vec <- msa_upd[,ps]
    obs <- unique(position_vec)
    
    matrix_type <- matrix(0, 5, 5)
    keep <- c()
    
    polarity1 <- length(intersect(obs, unlist(strsplit(type_aas[1,1], split = ","))))/length(unlist(strsplit(type_aas[1,1], split = ",")))
    polarity2 <- length(intersect(obs, unlist(strsplit(type_aas[1,2], split = ","))))/length(unlist(strsplit(type_aas[2,1], split = ",")))
    matrix_type[1,1:2] <- c(polarity1, polarity2)
    if (length(which(matrix_type[1,]>0))==1) {
      keep <- rbind(keep, c(1, which(matrix_type[1,]>0), max(matrix_type[1,])))
    }
    
    charge1 <- length(intersect(obs, unlist(strsplit(type_aas[2,1], split = ","))))/length(unlist(strsplit(type_aas[2,1], split = ",")))
    charge2 <- length(intersect(obs, unlist(strsplit(type_aas[2,2], split = ","))))/length(unlist(strsplit(type_aas[2,2], split = ",")))
    charge3 <- length(intersect(obs, unlist(strsplit(type_aas[2,3], split = ","))))/length(unlist(strsplit(type_aas[2,3], split = ",")))
    matrix_type[2,1:3] <- c(charge1, charge2, charge3)
    if (length(which(matrix_type[2,]>0))==1) {
      keep <- rbind(keep, c(2, which(matrix_type[2,]>0), max(matrix_type[2,])))
    }
    
    aromaticity1 <- length(intersect(obs, unlist(strsplit(type_aas[3,1], split = ","))))/length(unlist(strsplit(type_aas[3,1], split = ",")))
    aromaticity2 <- length(intersect(obs, unlist(strsplit(type_aas[3,2], split = ","))))/length(unlist(strsplit(type_aas[3,2], split = ",")))
    aromaticity3 <- length(intersect(obs, unlist(strsplit(type_aas[3,3], split = ","))))/length(unlist(strsplit(type_aas[3,3], split = ",")))
    matrix_type[3,1:3] <- c(aromaticity1, aromaticity2, aromaticity3)
    if (length(which(matrix_type[3,]>0))==1) {
      keep <- rbind(keep, c(3, which(matrix_type[3,]>0), max(matrix_type[3,])))
    }
    
    size1 <- length(intersect(obs, unlist(strsplit(type_aas[4,1], split = ","))))/length(unlist(strsplit(type_aas[4,1], split = ",")))
    size2 <- length(intersect(obs, unlist(strsplit(type_aas[4,2], split = ","))))/length(unlist(strsplit(type_aas[4,2], split = ",")))
    size3 <- length(intersect(obs, unlist(strsplit(type_aas[4,3], split = ","))))/length(unlist(strsplit(type_aas[4,3], split = ",")))
    matrix_type[4,1:3] <- c(size1, size2, size3)
    if (length(which(matrix_type[4,]>0))==1) {
      keep <- rbind(keep, c(4, which(matrix_type[4,]>0), max(matrix_type[4,])))
    }
    
    electronic1 <- length(intersect(obs, unlist(strsplit(type_aas[5,1], split = ","))))/length(unlist(strsplit(type_aas[5,1], split = ",")))
    electronic2 <- length(intersect(obs, unlist(strsplit(type_aas[5,2], split = ","))))/length(unlist(strsplit(type_aas[5,2], split = ",")))
    electronic3 <- length(intersect(obs, unlist(strsplit(type_aas[5,3], split = ","))))/length(unlist(strsplit(type_aas[5,3], split = ",")))
    electronic4 <- length(intersect(obs, unlist(strsplit(type_aas[5,4], split = ","))))/length(unlist(strsplit(type_aas[5,4], split = ",")))
    electronic5 <- length(intersect(obs, unlist(strsplit(type_aas[5,5], split = ","))))/length(unlist(strsplit(type_aas[5,5], split = ",")))
    matrix_type[5,1:5] <- c(electronic1, electronic2, electronic3, electronic4, electronic5)
    if (length(which(matrix_type[5,]>0))==1) {
      keep <- rbind(keep, c(5, which(matrix_type[5,]>0), max(matrix_type[5,])))
    }
    
    chosen_aas_aa <- c()
    if (length(keep)>3) {
      chosen <- keep[which(keep[,3]==max(keep[,3])),]
      if (length(chosen)>3){
        for (kl in 1:length(chosen[,1])) {
          chosen_aas_aa <- c(chosen_aas_aa, unlist(strsplit(type_aas[chosen[kl,1], chosen[kl,2]], split = ", ")))
        }
        chosen_aas_aa <- unique(chosen_aas_aa)
      } else {
        chosen_aas_aa <- unlist(strsplit(type_aas[chosen[1], chosen[2]], split = ", "))
      }
    } else if (length(keep)==3) {
      chosen <- keep
      chosen_aas_aa <- unlist(strsplit(type_aas[chosen[1], chosen[2]], split = ", "))
    } else {
      chosen <- 0
      chosen_aas_aa <- 0
    }
    
    diff_leaves <- matrix(0, num_leaves, 20)
    if (num_nodes == 1) {
      diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
      diff_nodes <- matrix(0, 1, 20)
      vect_human <- matrix_prob
    } else {
      diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
      diff_nodes <- matrix(0, num_nodes-1, 20)
      diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
      vect_human <- matrix_prob[node_human - num_leaves,]
    }
    diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
    
      score_aa <- matrix(0,1,20)
      
      if (num_nodes != 1) {
        s1 <- sapply(1:20, function(ii){
          if (num_nodes == 2) {
            dif_pr <- diff_nodes[ii]
          } else {
            dif_pr <- diff_nodes[1:(num_nodes-1),ii]
          }
          dif_pr[dif_pr<0] <- 0
          sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
          score_aa[ii] <<- score_aa[ii] + sum(weight_node[sel_node] * dif_pr)
        })
      }
      
      ### NOVEL 29.03
      aa_f <- position_num[human_plc]
      if (aa_f != 21) {
        vect_human[aa_f]<-0
      }
      
      score_aa_without_leaf <- score_aa
      score_aa_without_leaf <- score_aa_without_leaf + weight_node[(node_human-num_leaves)]*vect_human
      score_aa <- score_aa + weight_node[(node_human-num_leaves)]*vect_human
      
      s2 <- sapply(1:20, function(ii){
        diff_lf <- diff_leaves[1:num_leaves,ii]
        diff_lf[gaps] <-  0
        diff_lf[diff_lf<0] <- 0
        
        s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
        score_aa[ii] <<- score_aa[ii] + s1
      })
      
      aa_f <- position_num[human_plc]
      if (aa_f != 21){
        score_aa[aa_f] <- score_aa[aa_f] + weight_leaf[human_plc]*1
        score_aa_without_leaf[aa_f] <- score_aa_without_leaf[aa_f] + weight_leaf[human_plc]*1
      }
      
      sum_exc_max <- sum(score_aa_without_leaf)-max(score_aa_without_leaf)
      diversity <- (-(length(which(score_aa_without_leaf<0.0001))*0.1)/20+0.1)*(sum_exc_max)
      
      chosen_aas <- aa_to_num(chosen_aas_aa)
      if (length(chosen_aas)==1 &&  (chosen_aas==21 || chosen_aas_aa==0)) {
        score_aa_prop$score_aa <- score_aa/(num_nodes+num_leaves)
      } else {
        score_aa_prop$score_aa <- matrix(0, 1, 20)
        score_aa_prop$score_aa[chosen_aas] <- (score_aa[chosen_aas]*0.9 + diversity)/(num_nodes+num_leaves)
        score_aa_prop$score_aa[setdiff(1:20,chosen_aas)] <- (score_aa[-chosen_aas]*0.9)/(num_nodes+num_leaves)
      }
      score_aa_prop$div <- diversity
  }
  ############################## ML FEATURES ############################## 
  if (num_nodes == 1){
    human_node_prob <- matrix_prob
    prob_ordered <- matrix_prob
  } else {
    human_node_prob <- matrix_prob[(node_human - num_leaves),]
    prob_ordered <- matrix_prob[(as.numeric(names(sort(d_n))) - num_leaves), ]
  }
  
  if (length(matrix_prob)==20){
    subs_set <- names(matrix_prob)
    human_node_dist <- d_n
  } else {
    subs_set <- colnames(matrix_prob)
    human_node_dist <- d_n[as.character(node_human)]
  }
  
  acceptance_prob <- 0.1
  #at which nodes, a substitution is accepted (P(subs > acceptance prob))
  accepted_nodes <- list()
  binary_accepted <- list()
  n_of_nodes_accepted <- list()
  dist_accepted_min <- list()
  accepted_ratio <- list()
  first_observed_node <- list()
  S1 <- sapply(subs_set, function(aa){
    if (length(matrix_prob)==20){
      accepted_nodes[[aa]] <<- which(matrix_prob[aa] > acceptance_prob)
    } else {
      accepted_nodes[[aa]] <<- which(matrix_prob[,aa] > acceptance_prob)
    }
    binary_accepted[[aa]] <<- 1*(length(accepted_nodes[[aa]]) > 0)
    n_of_nodes_accepted[[aa]] <<- length(accepted_nodes[[aa]])
    dist_min <- min(d_n[accepted_nodes[[aa]]])
    dist_min <- ifelse(is.infinite(dist_min), 100, dist_min)
    dist_accepted_min[[aa]] <<- dist_min
    
    accepted_ratio[[aa]] <<- n_of_nodes_accepted[[aa]] / num_nodes
    sorted_dists <- sort(d_n)
    first_observed_node[[aa]] <<- min(which(dist_min == sorted_dists))
  })
  
  
  d_n_ordered <- sort(d_n)
  seqq <- ceiling(num_nodes * c(1:20)*0.05)
  
  
  if (num_nodes >= 2) {
    weighted_probs <- t(sapply(seqq, function(i){colSums(prob_ordered[(1:i), ] * matrix(exp(-d_n_ordered[1:i]), nrow = i, ncol = ncol(prob_ordered), byrow = FALSE)) / i}))
    rownames(weighted_probs) <- seqq
    colnames(weighted_probs) <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  } else {
    weighted_probs <- t(sapply(seqq, function(i){colSums(prob_ordered * matrix(exp(-d_n_ordered), nrow = i, ncol = 20, byrow = FALSE)) / i}))
    colnames(weighted_probs) <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  }
  
  if (ref_aa!="-"){
    if (num_nodes >= 2) {
      tree_kernel <- exp(-as.matrix(dist(t(matrix_prob))^2)/(2*(mean(as.matrix(dist(t(matrix_prob)))))))
      tree_kernel <- tree_kernel[ref_aa,]
    } else {
      tree_kernel <- matrix_prob
    }
  } else {
    tree_kernel <- matrix(0,1,20)
  }
  
  gap_freq <- sum(msa[,position] == "-") / (nrow(msa))
  
  if (num_leaves == 1) {
    pos_freq <- matrix(0,1,20)
    weighted_pos_freq <- matrix(0,1,20)
    if (aa_f != 21){
      pos_freq[aa_f] <- 1
      weighted_pos_freq[aa_f] <- 1
    }
    colnames(pos_freq) <- num_to_aa(1:20)
    colnames(weighted_pos_freq) <- num_to_aa(1:20)
  } else {
    pos_freq <- (sapply(c(num_to_aa(1:20)), function(a){sum(1*(msa_upd[,position] == a)) / sum(msa_upd[,position]!="-")}))
    weighted_pos_freq <- (sapply(c(num_to_aa(1:20)), function(a){sum(1*(msa_upd[(msa_upd[,position]!="-"),position] == a) * weight_leaf[(msa_upd[,position]!="-")])/sum(msa_upd[,position]!="-")}))
    
  }
  
  aa_class <- aaComp(ref_aa)
  aa_class <- aa_class[[1]][,1]
  
  ml_features <- list(ref_aa = ref_aa, binary_accepted = binary_accepted, n_of_nodes_accepted = n_of_nodes_accepted,
                      dist_accepted_min = dist_accepted_min, human_node_dist = human_node_dist, human_node_prob = human_node_prob,
                      weighted_probs = weighted_probs, tree_kernel = tree_kernel, aa_class = aa_class, pos_freq = pos_freq, gap_freq = gap_freq,
                      weighted_pos_freq = weighted_pos_freq, independent_increments_leaves = ind_inc_leaves, independent_increments_nodes = ind_inc_nodes,
                      independent_increments_total = ind_inc_total, accepted_ratio = accepted_ratio, first_observed_node = first_observed_node,
                      score_aaprop = score_aa_prop$score_aa, score_aadiv = score_aa_prop$div)

  #################################
  ########AllDifferences with CN3########
  diff_leaves <- matrix(0, num_leaves, 20)
  if (num_nodes == 1) {
    diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
    diff_nodes <- matrix(0, 1, 20)
    vect_human <- matrix_prob
  } else {
    diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
    diff_nodes <- matrix(0, num_nodes-1, 20)
    diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
    vect_human <- matrix_prob[node_human - num_leaves,]
  }
  diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
  
  score_alldiff <- matrix(0,1,20)
  if (num_nodes != 1) {
    s1 <- sapply(1:20, function(ii){
      if (num_nodes == 2) {
        dif_pr <- diff_nodes[ii]
      } else {
        dif_pr <- diff_nodes[1:(num_nodes-1),ii]
      }
      sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
      score_alldiff[ii] <<- score_alldiff[ii] + sum(weight_node[sel_node] * dif_pr)
    })
  }

  score_alldiff_without_leaf <- score_alldiff
  score_alldiff_without_leaf <- score_alldiff_without_leaf + weight_node[(node_human-num_leaves)]*vect_human
  score_alldiff <- score_alldiff + weight_node[(node_human-num_leaves)]*vect_human
  
  s2 <- sapply(1:20, function(ii){
    diff_lf <- diff_leaves[1:num_leaves,ii]
    diff_lf[gaps] <-  0
    
    s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
    score_alldiff[ii] <<- score_alldiff[ii] + s1
  })
  
  aa_f <- position_num[human_plc]
  if (aa_f != 21){
    score_alldiff[aa_f] <- score_alldiff[aa_f] + weight_leaf[human_plc]*1
    score_alldiff_without_leaf[aa_f] <- score_alldiff_without_leaf[aa_f] + weight_leaf[human_plc]*1
  }
  
  colnames(score_alldiff) <- num_to_aa(1:20)
  score_alldiff_norm <- score_alldiff/(num_nodes+num_leaves)
  
  ml_features$allDiff_wl <- score_alldiff_norm
  
  ########AllProbwith CN3########
  score_allprob <- matrix(0,1,20)
  
  if (num_leaves == 1) {
    if (aa_f != 21){
      score_allprob[aa_f] <- 1
    }
    score_allprob <- score_allprob + matrix_prob*weight_node
    
  } else {
    diff_leaves <- matrix(0, num_leaves, 20)
    diff_leaves <- prob_leaves  
    
    if (num_nodes == 1) {
      diff_nodes <- matrix_prob
    } else {
      diff_nodes <- matrix(0, num_nodes-1, 20)
      diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ]
    }
    
    s1 <- sapply(1:20, function(ii){
      if (length(diff_nodes)==20) {
        dif_pr <- diff_nodes[ii]
        score_allprob[ii] <<- score_allprob[ii] + sum(weight_node * dif_pr)
      } else {
        dif_pr <- diff_nodes[1:(num_nodes-1),ii]
        sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
        score_allprob[ii] <<- score_allprob[ii] + sum(weight_node[sel_node] * dif_pr)
      }
    })
    
    s2 <- sapply(1:20, function(ii){
      diff_lf <- diff_leaves[1:num_leaves,ii]
      s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
      score_allprob[ii] <<- score_allprob[ii] + s1
    })
    
    if (aa_f != 21){
      score_allprob[aa_f] <- score_allprob[aa_f] + weight_leaf[human_plc]*1
    }
  }
  
  colnames(score_allprob) <- num_to_aa(1:20)
  score_allprob_norm <- score_allprob/(num_nodes+num_leaves)
  
  ml_features$allProb_wl <- score_allprob_norm
  
  #################################
  #################################
  ############################## END ############################## 
  ml_features_prot_scale <- list()
  ml_features_prot_scale$ml_features <- ml_features
  ml_features_prot_scale$protscale_scores <- scale_scores
  ml_features_prot_scale$protscale_scores_wo_phact <- scale_scores_wo_phact
  
  
  return(ml_features_prot_scale)
}


