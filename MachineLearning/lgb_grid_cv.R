lgb_grid_cv <- function(X, y, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set) {
  nrow_mat <- length(eta_set)*length(max_depth_set)*length(num_leaves_set[[1]])*length(min_data_set)*length(lambda_l1_set)*length(lambda_l2_set)*length(min_gain_to_split_set)*length(bagging_fraction_set)*length(bagging_freq_set)*length(feature_fraction_set)*length(min_hessian_set)*length(path_smooth_set)
  auroc_matrix <- array(NA, dim = c(nrow_mat,15,fold_count), dimnames = list(c(1:nrow_mat), c("train_auc", "auc", "eta", "depth", "n_leaves", "min_data", "lambda_l1", "lambda_l2", "min_gain", "bagging_frac", "bagging_freq", "feature_frac", "min_hessian", "path_smooth", "best_nrounds"), c(1:fold_count)))
  
  for (fold in 1:fold_count) {
    train_indices <- fold_indices_train[[fold]]
    test_indices <- fold_indices_test[[fold]]
    X_train <- as.matrix(X[train_indices,])
    X_test <- as.matrix(X[test_indices,])
    X_train <- scale(X_train)
    X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
    
    y_train <- y[train_indices]
    y_test <- y[test_indices]
    
    dtrain <- lgb.Dataset(data = as.matrix(X_train), label = 1*(y_train==1))
    dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = as.matrix(X_test), label = 1*(y_test==1))
    
    print(sprintf("running replication = %d, fold = %d", replication, fold))
    i <- 0
    for (eta in eta_set) {
      for (max_depth in max_depth_set) {
        for (num_leaves in num_leaves_set[[sprintf("%s",max_depth)]]) {
          for (min_data in min_data_set) {
            for (lambda_l1 in lambda_l1_set) {
              for (lambda_l2 in lambda_l2_set) {
                for (min_gain_to_split in min_gain_to_split_set) {
                  for (bagging_fraction in bagging_fraction_set) {
                    for (bagging_freq in bagging_freq_set) {
                      for (feature_fraction in feature_fraction_set) {
                        for (min_hessian in min_hessian_set) {
                          for (path_smooth in path_smooth_set) {
                            t1 <- Sys.time()
                            i <- i+1
                            parameters <- list()
                            parameters$learning_rate <- eta
                            parameters$max_depth <- max_depth
                            parameters$num_leaves <- num_leaves
                            parameters$min_data <- min_data
                            parameters$lambda_l1 <- lambda_l1
                            parameters$lambda_l2 <- lambda_l2
                            parameters$min_gain_to_split <- min_gain_to_split
                            parameters$bagging_fraction <- bagging_fraction
                            parameters$bagging_freq <- bagging_freq
                            parameters$feature_fraction <- feature_fraction
                            parameters$min_hessian <- min_hessian
                            parameters$path_smooth <- path_smooth
                            parameters$boosting <- boosting_type
                            # parameters$obj <- as.vector("binary")
                            parameters$metric <- as.vector("auc")
                            # parameters$is.unbalance <-  T
                            # parameters$histogram_pool_size <- 8*1024
                            parameters$scale_pos_weight <- sum(y_train == -1) / sum(y_train == 1)
                            parameters$verbosity <-  -1
                            state <- lgb.train(params = parameters, data = dtrain, obj = "binary", nrounds = nrounds, valids = list(train = dtrain, test = dtest),
                                               early_stopping_rounds = ceiling(nrounds * 0.05), set.seed(1), eval_freq = 50, reset_data = T)          
                            gc(verbose = FALSE)
                            # prediction <- predict(state, as.matrix(X_test))
                            # train_pred <- predict(state, as.matrix(X_train))
                            # train_auroc <- auc(roc(train_pred, as.factor(1 * (y_train == +1))))
                            # auroc <- auc(roc(prediction, as.factor(1 * (y_test == +1))))
                            train_auroc <- state$record_evals$train$auc$eval[[state$best_iter]]
                            test_auroc <- state$record_evals$test$auc$eval[[state$best_iter]]
                            auroc_matrix[i, ,fold] <- c(train_auroc, test_auroc, eta, max_depth, num_leaves, min_data, 
                                                        lambda_l1, lambda_l2, min_gain_to_split, bagging_fraction, bagging_freq, 
                                                        feature_fraction, min_hessian, path_smooth, state$best_iter)
                            t2 <- Sys.time()
                            print(t2-t1)
                            
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(auroc_matrix)
}
