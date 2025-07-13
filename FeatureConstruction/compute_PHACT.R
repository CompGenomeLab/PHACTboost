position_score <- function(ps, x, msa, human_id, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes_prev, nodes_raxml_prev, num_leaves_prev, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, parameters) {
  position <- ps
  #print(ps)  
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
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  
  trims <- unique(trims)
  if (length(grep(human_id, trims))>0){
    ALL_SCORES <- list()
    pl <- 1
    for (parameter in parameters) {
      scores <- list()
      scores$score_with_leaf <- matrix(0,1,20)
      scores$score_without_leaf <- matrix(0,1,20)
      scores$diversity <- 0
      scores$score_nogap <- matrix(0,1,20)
      
      ALL_SCORES <- c(ALL_SCORES, scores)
      names(ALL_SCORES)[pl:(pl+3)] <- paste(parameter, "_", names(ALL_SCORES)[pl:(pl+3)], sep = "")
      pl <- pl + 4
    }
    
    ALL_SCORES$trims <- num_leaves
    ALL_SCORES$div_pos <- 0
    ALL_SCORES$div_pos2 <- 0
  } else {
    if (length(trims)>=1) {
      tree_new <- drop.tip(tr_org, trims)
      tree_new_info <- as.data.frame(as_tibble(tree_new))
      
      num_leaves <- length(tree_new[["tip.label"]])
      num_nodes <- tree_new[["Nnode"]]
      
      nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_new_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
      names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
      
      elim_nodes <- setdiff(nodes_raxml_prev, nodes_raxml)
      if (length(elim_nodes)>0){
        els <- t(mapply(function(jj){which(as.numeric(nodes_raxml_prev) == elim_nodes[jj])}, rep(1:length(elim_nodes))))
        matrix_prob <- matrix_prob[-els,]
      }
      
      if (length(trims)>=1){
        msa_upd <- msa_upd[-t(mapply(function(jj){which(rownames(msa_upd)==trims[jj])}, rep(1:length(trims)))),]
      }
    } else {
      tree_new <- tr_org
      tree_new_info <-  tree_info
    }
    
    ALL_SCORES <- list()
    if (num_leaves == 1){
      score <- matrix(0,1,20)
      aa_f <- aa_to_num(msa_upd[ps])
      if (aa_f!=21){
        score[aa_f] <- 1
      }
      diversity <- 0
      position_vec <- msa_upd[ps]
      position_num <- aa_to_num(position_vec)
      
      pl <- 1
      for (parameter in parameters){
        scores <- list()
        scores$score_with_leaf <- score/(num_nodes+num_leaves)
        scores$score_without_leaf <- score/num_nodes
        
        scores$diversity <- diversity
        num_gaps <- sum(position_num==21)
        scores$score_nogap <- score/(num_nodes+num_leaves-num_gaps)
        
        ALL_SCORES <- c(ALL_SCORES, scores)
        names(ALL_SCORES)[pl:(pl+3)] <- paste(parameter, "_", names(ALL_SCORES)[pl:(pl+3)], sep = "")
        pl <- pl + 4
      }
      
    } else {
      h_name <- human_id
      human_codeml <- names_all[grep(pattern = h_name, x = names_all, fixed = TRUE)]
      leaf_human <- tree_new_info[which(tree_new_info$label == human_codeml), "node"]
      human_plc <- leaf_human
      node_human <- tree_new_info[which(tree_new_info$label == human_codeml), "parent"]
      
      
      dd_node <- dist.nodes(tree_new)
      dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
      dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
      
      # Connections between leaves & nodes
      chosen_leaves <- tree_new_info[1:num_leaves,c("parent", "node")]
      # Connections between nodes & nodes
      chosen_nodes <- tree_new_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
      leaf_names <- tree_new_info$label
      
      human_leaf_len <- as.double(tree_new_info[human_plc, "branch.length"])
      
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
      
      
      position_vec <- msa_upd[, ps]
      
      position_num <- aa_to_num(position_vec)
      prob_leaves <- matrix(0, num_leaves, 20)
      prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
      
      gaps <- which(position_num == 21)
      
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
      
      pl <- 1
      ################## weights
      for (parameter in parameters) {
        weights <- weight_fnc(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn)
        weight_leaf <- weights[1:num_leaves]
        weight_node <- tail(weights,num_nodes)
        
        score <- matrix(0,1,20)
        
        if (num_nodes != 1) {
          s1 <- sapply(1:20, function(ii){
            if (num_nodes == 2) {
              dif_pr <- diff_nodes[ii]
            } else {
              dif_pr <- diff_nodes[1:(num_nodes-1),ii]
            }
            dif_pr[dif_pr<0] <- 0
            sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
            score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr)
          })
        }
        
        ### NOVEL 29.03
        aa_f <- position_num[human_plc]
        if (aa_f != 21) {
          vect_human[aa_f]<-0
        }
        
        score_without_leaf <- score
        score_without_leaf <- score_without_leaf + weight_node[(node_human-num_leaves)]*vect_human
        score <- score + weight_node[(node_human-num_leaves)]*vect_human
        
        s2 <- sapply(1:20, function(ii){
          diff_lf <- diff_leaves[1:num_leaves,ii]
          diff_lf[gaps] <-  0
          diff_lf[diff_lf<0] <- 0
          
          s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
          score[ii] <<- score[ii] + s1
        })
        
        aa_f <- position_num[human_plc]
        if (aa_f != 21){
          score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*1
          score_without_leaf[aa_f] <- score_without_leaf[aa_f] + weight_leaf[human_plc]*1
        }
        
        sum_exc_max <- sum(score_without_leaf)-max(score_without_leaf)
        diversity <- (-(length(which(score_without_leaf<0.0001))*0.1)/20+0.1)*(sum_exc_max)
        
        scores <- list()
        scores$score_with_leaf <- (score*0.9 + diversity)/(num_nodes+num_leaves)
        scores$score_without_leaf <- (score_without_leaf*0.9 + diversity)/num_nodes
        
        scores$diversity <- diversity
        num_gaps <- sum(position_num==21)
        scores$score_nogap <- (score*0.9 + diversity)/(num_nodes+num_leaves-num_gaps)
        
        ALL_SCORES <- c(ALL_SCORES, scores)
        names(ALL_SCORES)[pl:(pl+3)] <- paste(parameter, "_", names(ALL_SCORES)[pl:(pl+3)], sep = "")
        pl <- pl + 4
      }
      
    }
    
    ALL_SCORES$trims <- length(trims)/num_leaves
    ALL_SCORES$div_pos <- length(unique(position_vec))
    if (ps==1) {
      un1 <- length(unique(msa[which(msa[,(ps+1)]!="-"),(ps+1)])) + length(unique(position_vec))
    } else if (ps == length(msa[1,])) {
      un1 <- length(unique(msa[which(msa[,(ps-1)]!="-"),(ps-1)])) + length(unique(position_vec))
    } else {
      un1 <- length(unique(msa[which(msa[,(ps-1)]!="-"),(ps-1)])) + length(unique(position_vec)) +
        length(unique(msa[which(msa[,(ps+1)]!="-"),(ps+1)]))
    }
    ALL_SCORES$div_pos2 <- un1
    
  }
  
  
  return(ALL_SCORES)
}
