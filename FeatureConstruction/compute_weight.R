weight_fnc <- function(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn, mxx) {
  # print(parameter)
  if (parameter=="0"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
  } else if (parameter=="0_MinNode"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
  } else if (parameter=="0_MinNode_Mix"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node1 <- 1/d_n2
    weight_leaf1 <- 1/d_l2
    
    param <- mean(c(d_n, d_l))
    weight_leaf2 <- exp(-d_l^2/param^2)
    weight_node2 <- exp(-d_n^2/param^2)
    
    weight_node<-sqrt(weight_node1*weight_node2)
    weight_leaf<-sqrt(weight_leaf1*weight_leaf2)
    
  } else if (parameter=="0_MinNode_Mix2"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node1 <- 1/d_n2
    weight_leaf1 <- 1/d_l2
    
    param <- mean(c(d_n, d_l))
    weight_leaf2 <- exp(-d_l^2/param^2)/2
    weight_node2 <- exp(-d_n^2/param^2)/2
    
    weight_node<-sqrt(weight_node1*weight_node2)
    weight_leaf<-sqrt(weight_leaf1*weight_leaf2)
    
  } else if (parameter == "mean"){
    param <- mean(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  } else if (parameter == "median"){
    param <- median(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  } else if (parameter == "X"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    min_n <- min(d_n)
    min_ch <- min(min_l, min_n)
    d_l2 <- d_l-min_ch
    d_n2 <- d_n-min_ch
    weight_leaf <- exp(-d_l^2)/2
    weight_node <- exp(-d_n^2)/2
    weight_leaf[human_plc] <- 1
  } else if (parameter == "CountNodes_1"){      ### NewFunction 1 (15Nov)
    weight_node = (exp(-d_n^2) + 1/nodes_conn)/2
    weight_leaf = (exp(-d_l^2) + 1/leaves_conn)/2
  } else if (parameter == "CountNodes_2"){      ### NewFunction 2 (15Nov)
    weight_node = (exp(-d_n^2) + exp(-nodes_conn^2))/2
    weight_leaf = (exp(-d_l^2) + exp(-leaves_conn^2))/2
  } else if (parameter == "CountNodes_3"){      ### NewFunction 3 (15Nov)
    weight_node = sqrt(exp(-d_n^2)*1/nodes_conn)
    weight_leaf = sqrt(exp(-d_l^2)*1/leaves_conn)
  } else if (parameter == "CountNodes_4"){      ### NewFunction 4 (15Nov)
    weight_node = exp(-(sqrt(d_n*nodes_conn))^2)
    weight_leaf = exp(-(sqrt(d_l*leaves_conn))^2)
  } else if (parameter == "CountNodes_5"){
    param <- mean(c(d_n, d_l))
    weight_node = sqrt(exp(-d_n^2/param^2)*1/nodes_conn)
    weight_leaf = sqrt(exp(-d_l^2/param^2)*1/leaves_conn)
  } else if (parameter == "Equal"){
    
    weight_leaf <- matrix(1,1,length(d_l))
    weight_node <- matrix(1,1,length(d_n))
    
  } else if (parameter == "MinThreshold"){
    max_dis <- max(c(d_l, d_n))
    param_min <- as.double(param_min)    
    weight_node <- (-1+param_min)/max_dis*d_n + 1
    weight_leaf <- (-1+param_min)/max_dis*d_l + 1
  } else if (parameter == "MinThreshold_Gauss"){
    max_dis <- max(c(d_l, d_n))
    param_min <- as.double(param_min)
    param <- sqrt(-max_dis^2/log(param_min))
    
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
    
  } else {
    param <- as.double(parameter)
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  }
  weights = c(weight_leaf, weight_node)
  
  
  return(weights)
}
