rm(list=ls())

library(lightgbm)
library(AUC)

ids <- args[1]
codon_info_path <- args[2]
train_path <- args[3]
input_features_path <- args[4]
final_model_path <- args[5]

for (index in 1:length(ids)){
  id <- ids[index]

  # Load training set
  load(train_path)

  # Load input features data
  load(sprintf("%s/%s.RData", input_features_path, id))
  
  add_pos <- matrix(0, length(data$UNIPROTKB), 1)
  add_snp <- matrix(0, length(data$UNIPROTKB), 1)
  
  add_pos[which(data$Positions==1)] <- 1
  
  codon_info <- readxl::read_xlsx(codon_info_path)
  
  for (i in 1:length(data$Alt_AA)){
    snp <- 0
    r <- data$Ref_AA[i]
    l <- which(codon_info$Column2==r)
    a <- data$Alt_AA[i]
    if (length(grep(a, codon_info$Column3[l]))>0){
      snp <- 1
    } else if (length(grep(a, codon_info$Column4[l]))>0) {
      # print(i)
      snp <- 2
    }
    add_snp[i] <- snp
  }
  
  nSNP <- add_snp
  data <- cbind(data, nSNP)
  data <- cbind(data, add_pos)
  colnames(data)[492] <- "Pos1"
  
  replication <- 1
  phact_choice <- "with_phact"
  param_choice <- "CountNodes_3"
  star_choice <- "with_0"
  eliminate_circularity <- ""
  
  set.seed(1903 * replication)
  
  all_train <- train
  all_test <- data
  
  y_train <- all_train$variant_info
  xx <- which(colnames(all_train)=="vars")
  X_train <- all_train[, -c(1:11, xx, grep("SIFT", colnames(all_train)))]
  X_train <- X_train[, (!grepl("wol", colnames(X_train)))]
  
  y_test <- all_test$variant_info
  xx2 <- which(colnames(all_test)=="vars")
  X_test <- all_test[, -c(1:4, xx2, grep("SIFT", colnames(all_test)))]
  X_test <- X_test[, (!grepl("wol", colnames(X_test)))]
  
  elims <- which(apply(X_train, 2, sd) == 0)
  X_train <- X_train[,-elims]
  X_test <- X_test[,-elims]
  
  parameters <- c("0.1", "0.5", "1", "2", "3",  "5", "mean", "median", "0", paste0("CountNodes_", c(1:5)), "0_MinNode_Mix", "0_MinNode_Mix2", "Equal", "max05_Gauss")
  if (param_choice  %in% parameters) {
    left_params <- parameters[parameters != param_choice]
    param_cols <- which(grepl("param", colnames(X_train)))
    param_selected <- param_cols[which(!grepl(paste0("param_", left_params, collapse = "|"), colnames(X_train)[param_cols]))]
    
    selected_cols <- c(param_selected, which(!grepl("param", colnames(X_train))))
    X_train <- X_train[, selected_cols]
    X_test <- X_test[, selected_cols]
  }
  
  X_train <- scale(X_train)
  
  train_center <- attr(X_train, "scaled:center")
  train_scale <- attr(X_train, "scaled:scale")
  
  train_info <- list(center = train_center, scale = train_scale, features = colnames(X_train))
  save(train_info, file = "train_info.RData")
  
  X_test <- X_test[,train_info$features]
  X_test <- (X_test - matrix(train_info$center, nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(train_info$scale, nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
  colnames(X_test) <- train_info$features
  
  state <- lgb.load(final_model_path)
  PhactBoost_Scores <- predict(state, as.matrix(X_test), params = list(predict_disable_shape_check=T))
  data <- cbind(PhactBoost_Scores, data)

  save(data, file = sprintf("PHACTboost_%s.RData", id))
}

