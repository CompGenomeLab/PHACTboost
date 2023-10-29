library(lightgbm)
library(AUC)
source("./lgb_grid_cv_18082022.R")

load("TrainSet.RData")
load("TestSet.RData")

all_train <- train
all_test <- test

args <- commandArgs(trailingOnly = TRUE)
replication <- args[[1]]

param_choice <- args[[2]]
rpath <- args[[3]]
result_path <- sprintf("./%s_%s_%s_%s", rpath, param_choice, replication)

if(dir.exists(sprintf("%s", result_path)) == FALSE) {
  dir.create(sprintf("%s",result_path))
}

set.seed(1903 * replication)

y_train <- all_train$variant_info
el1 <- which(colnames(all_train)=="vars")
X_train <- all_train[, -c(1:11, el1, grep("SIFT", colnames(all_train)))]
X_train <- X_train[, (!grepl("wol", colnames(X_train)))]
elims <- which(apply(X_train, 2, sd) == 0)
X_train <- X_train[,-elims]

y_test <- all_test$variant_info
el2 <- which(colnames(all_test)=="vars")
X_test <- all_test[, -c(1:11, el2, grep("SIFT", colnames(all_test)))]
X_test <- X_test[, (!grepl("wol", colnames(X_test)))]
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

epsilon <- 1e-3
fold_count <- 4
train_ratio <- 0.8
boosting_type <- "gbdt"

train_negative_indices <- which(y_train == -1)
train_positive_indices <- which(y_train == 1)

positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))
negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))

fold_indices_train <- list()
fold_indices_test <- list()
for (fold in 1:fold_count) {
  train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
  test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
  
  fold_indices_train[[fold]] <- train_indices
  fold_indices_test[[fold]] <- test_indices
}

features <- colnames(X_train)
features <- gsub(",", "_", features)
colnames(X_train) <- features
colnames(X_test) <- features

phact <- 1
wp <- 1
ref <- 1
delta <- 1
aa_class <- 1 
expasy_phact <- 1 
expasy_wophact <- 1
expasy_only  <- 1 
msa <- 1   
phytree <- 1

feature_groups <- list(phact = phact, wp = wp, ref = ref, delta = delta, aa_class = aa_class, expasy_phact = expasy_phact,
                       expasy_wophact = expasy_wophact, expasy_only = expasy_only, msa = msa, phytree = phytree)

selected_features <- colnames(X_train)

nrounds <- 5000
eta_set <- c(seq(0.001,0.009,0.001), seq(0.01,0.3,0.01))
max_depth_set <- 6
num_leaves_set <- list()
num_leaves_set[[sprintf("%s", max_depth_set)]] <- 50
min_data_set <- 20
lambda_l1_set <- 0
lambda_l2_set <- 0
min_gain_to_split_set <- 0
bagging_fraction_set <- 0.7
bagging_freq_set <- 0
feature_fraction_set <- 0.7
min_hessian_set <- 1e-3
path_smooth_set <- 0

print("Tuning eta")
auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

eta_star <- best_auroc[["eta"]]
if(eta_star >= 0.01){
  eta_set <- c(eta_star - 0.005, eta_star, eta_star + 0.005)
}else{
  eta_set <- c(eta_star - 0.0005, eta_star, eta_star + 0.0005)
} 

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

eta_star <- best_auroc[["eta"]]
eta_set <- eta_star
print("Tuning max depth and num_leaves")
max_depth_set <- seq(3,7,2)
num_leaves_set <- list()
for (d in max_depth_set){
  leaf_set <- ceiling((2^d) * c(0.6,0.8,1))
  num_leaves_set[[sprintf("%s",d)]] <- leaf_set
}

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

max_depth_star <- best_auroc["depth"]
max_depth_set <- c((max_depth_star - 1): (max_depth_star + 1)) 

num_leaves_set <- list()
for (d in max_depth_set){
  leaf_set <- ceiling((2^d) * c(0.6,0.8,1))
  num_leaves_set[[sprintf("%s",d)]] <- leaf_set
}

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

max_depth_star <- best_auroc["depth"]
num_leaves_star <- best_auroc["n_leaves"]
max_depth_set <- max_depth_star
coef <- c(0.6,0.8,1)
coef_star <- coef[which(ceiling((2^max_depth_star)*coef) == num_leaves_star)]
num_leaves_set <- list()
num_leaves_set[[sprintf("%s", max_depth_star)]] <- unique(ceiling((2^max_depth_star)*c(coef_star-0.1, coef_star, coef_star+0.1)))
auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

num_leaves_star <- best_auroc["n_leaves"]
num_leaves_set <- list()
num_leaves_set[[sprintf("%s", max_depth_star)]] <- num_leaves_star

print("Tuning min data per leaf")
min_data_set <- c(seq(10,150,10))

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

min_data_star <-  best_auroc["min_data"]
min_data_set <- seq(min_data_star - 5, min_data_star + 5, 5)
min_data_set <- min_data_set[min_data_set > 0]

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

min_data_star <-  best_auroc["min_data"]
min_data_set <- min_data_star

print("Tuning min_hessian")
min_hessian_set <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

min_hessian_star <- best_auroc["min_hessian"]
min_hessian_set <- c(min_hessian_star*0.5, min_hessian_star, min_hessian_star*5)

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

min_hessian_star <- best_auroc["min_hessian"]
min_hessian_set <- min_hessian_star

print("Tuning path_smooth")

path_smooth_set <- c(seq(0,10,1))

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

path_smooth_star <- best_auroc["path_smooth"]
path_smooth_set <- c(path_smooth_star - 0.5, path_smooth_star - 0.25, path_smooth_star, path_smooth_star + 0.25, path_smooth_star + 0.5)
path_smooth_set <- path_smooth_set[path_smooth_set >= 0]

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

path_smooth_star <- best_auroc["path_smooth"]
path_smooth_set <- path_smooth_star

print("Tuning l1 regularization")

lambda_l1_set <- c(seq(0,20,5))

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

lambda_l1_star <- best_auroc["lambda_l1"]
if(lambda_l1_star == 0){
  lambda_l1_set <- c(seq(0,1,0.2), 2.5)
}else{
  lambda_l1_set <- seq(lambda_l1_star - 2.5, lambda_l1_star + 2.5, 2.5)
  lambda_l1_set <- lambda_l1_set[lambda_l1_set >= 0]
}

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

lambda_l1_star <- best_auroc["lambda_l1"]
lambda_l1_set <- lambda_l1_star

print("Tuning l2 regularization")

lambda_l2_set <- c(seq(0,20,5))

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

lambda_l2_star <- best_auroc["lambda_l2"]
if(lambda_l2_star == 0){
  lambda_l2_set <- c(seq(0,1,0.2), 2.5)
}else{
  lambda_l2_set <- seq(lambda_l2_star - 2.5, lambda_l2_star + 2.5, 2.5)
  lambda_l2_set <- lambda_l2_set[lambda_l2_set >= 0]
}

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

lambda_l2_star <- best_auroc["lambda_l2"]
lambda_l2_set <- lambda_l2_star

print("Tuning min gain to split")
min_gain_to_split_set <- c(seq(0,20,5))
auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)

mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

min_gain_to_split_star <- best_auroc["min_gain"]
if(min_gain_to_split_star == 0){
  min_gain_to_split_set <- c(seq(0,1,0.2), 2.5)
}else{
  min_gain_to_split_set <- seq(min_gain_to_split_star - 2.5, min_gain_to_split_star + 2.5, 2.5)
  min_gain_to_split_set <- min_gain_to_split_set[min_gain_to_split_set >= 0]
}

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

min_gain_to_split_star <- best_auroc["min_gain"]
min_gain_to_split_set <- min_gain_to_split_star

print("Tuning bagging fractions and bagging frequency")
bagging_fraction_set <- c(seq(0.5,1,0.1))
bagging_freq_set <- c(25, 50, 75, 100)

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

bagging_fraction_star <- best_auroc["bagging_frac"]
bagging_freq_star <- best_auroc["bagging_freq"]
bagging_freq_set <- bagging_freq_star

bagging_fraction_set <- c(bagging_fraction_star - 0.05, bagging_fraction_star, bagging_fraction_star + 0.05)
bagging_fraction_set <- bagging_fraction_set[bagging_fraction_set <= 1]

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

bagging_fraction_star <- best_auroc["bagging_frac"]
bagging_fraction_set <- bagging_fraction_star

print("Tuning feature fractions")
feature_fraction_set <- c(seq(0.1,1,0.1))

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

feature_fraction_star <- best_auroc["feature_frac"]
feature_fraction_set <- c(feature_fraction_star - 0.05, feature_fraction_star, feature_fraction_star + 0.05)
feature_fraction_set <- feature_fraction_set[feature_fraction_set <= 1]

auroc_matrix <- lgb_grid_cv(X_train, y_train, fold_count, fold_indices_train, fold_indices_test, boosting_type, eta_set, nrounds, max_depth_set, num_leaves_set, min_data_set , lambda_l1_set, lambda_l2_set, min_gain_to_split_set, bagging_fraction_set, bagging_freq_set, feature_fraction_set, min_hessian_set, path_smooth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]

feature_fraction_star <- best_auroc["feature_frac"]
feature_fraction_set <- feature_fraction_star

nrounds_star <- nrounds

print("Training the best model")

parameters <- list()
parameters$learning_rate <- eta_star
parameters$max_depth <- max_depth_star
parameters$num_leaves <- num_leaves_star
parameters$min_data <- min_data_star
parameters$lambda_l1 <- lambda_l1_star
parameters$lambda_l2 <- lambda_l2_star
parameters$min_gain_to_split <- min_gain_to_split_star
parameters$bagging_fraction <- bagging_fraction_star
parameters$bagging_freq <- bagging_freq_star
parameters$feature_fraction <- feature_fraction_star
parameters$min_hessian <- min_hessian_star
parameters$path_smooth <- path_smooth_star
parameters$boosting <- boosting_type
parameters$metric <- as.vector("auc")
parameters$scale_pos_weight <- sum(y_train == -1) / sum(y_train == 1)
parameters$verbosity <-  -1

X_train <- scale(X_train)
X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
colnames(X_test) <- colnames(X_train)
dtrain <- lgb.Dataset(data = as.matrix(X_train), label = 1*(y_train==1))
dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = as.matrix(X_test), label = 1*(y_test==1))

state <- lgb.train(params = parameters, data = dtrain, obj = "binary", nrounds = nrounds_star, valids = list(train = dtrain, test = dtest),
                   early_stopping_rounds = ceiling(nrounds_star * 0.05), set.seed(1), eval_freq = 50, reset_data = T)          
lgb.save(state, sprintf("%s/lightgbm_replication_%d_model.txt", result_path, replication))

train_auroc1 <- state$record_evals$train$auc$eval[[state$best_iter]]
test_auroc1 <- state$record_evals$test$auc$eval[[state$best_iter]]

train_prediction <- predict(state, as.matrix(X_train))
train_auroc <- auc(roc(train_prediction, as.factor(1 * (y_train == 1))))

test_prediction <- predict(state, as.matrix(X_test))
test_auroc <- auc(roc(test_prediction, as.factor(1 * (y_test == 1))))

prediction <- list(train_prediction = train_prediction, test_prediction = test_prediction, parameters = parameters, nrounds = nrounds_star,
                   selected_features = selected_features, feature_groups = feature_groups)

result <- list(train_auroc = train_auroc, test_auroc = test_auroc)
save("prediction", file = sprintf("%s/lightgbm_replication_%d_prediction.RData", result_path, replication))
save("result", file = sprintf("%s/lightgbm_replication_%d_result.RData", result_path, replication))

