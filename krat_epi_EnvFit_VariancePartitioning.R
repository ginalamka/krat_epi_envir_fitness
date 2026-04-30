# temp ridge regressions for variance partitioning - krat env
# # code from: https://chatgpt.com/share/68f00e1b-693c-8012-a353-89c657f46360

# CONSIDER IF INBREEDING/mkins_all VALUE SHOULD BE ADDED TO GENETIC (note many 0s for inbreeding)

## CODE FOR LRS
# Data subsets
methyl_only <- data.frame(methyl_subset, off_survive = sig_methyl_data$off_survive)
gen_only <- data.frame(sig_methyl_data[, c("lineage", "pop", "off_survive")])
env_only <- data.frame(sig_methyl_data[, c("month", "year", "rain", "temp", "totpopsz", "off_survive")])

# Combined (already defined)
combine_data_RRS <- data.frame(methyl_subset, sig_methyl_data[,c("lineage", "pop", "totpopsz",
                                                                 "month", "year", "rain", "temp",
                                                                 "off_survive")])

winds <- data[,35:ncol(data)]
all_data_RRS <- data.frame(winds, meta[,c("lineage", "pop", "totpopsz",
                                                             "month", "year", "rain", "temp",
                                                             "off_survive")])

library(caret)
library(ridge)

ridge_train <- function(df) {
  caret::train(off_survive ~ .,
               data = df,
               method = "glmnet",
               tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(2, 3, length = 100)),
               trControl = trainControl(method = "repeatedcv", number = 10, repeats = 100))
}

RRS_14k <- ridge_train(all_data_RRS)
RRS_all  <- ridge_train(combine_data_RRS)
RRS_meth <- ridge_train(methyl_only)

ridge_train <- function(df) {
  caret::train(off_survive ~ .,
               data = df,
               method = "glmnet",
               tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(0, 1, length = 100)),
               trControl = trainControl(method = "repeatedcv", number = 10, repeats = 100))
}
  
RRS_env  <- ridge_train(env_only)

ridge_train <- function(df) {
  caret::train(off_survive ~ .,
               data = df,
               method = "glmnet",
               tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(1, 2, length = 100)),
               trControl = trainControl(method = "repeatedcv", number = 10, repeats = 100))
}
RRS_gen  <- ridge_train(gen_only)

lambda_14k  <- RRS_14k$bestTune$lambda # 335
lambda_all  <- RRS_all$bestTune$lambda # 259
lambda_meth <- RRS_meth$bestTune$lambda # 319
lambda_gen  <- RRS_gen$bestTune$lambda # 32
lambda_env  <- RRS_env$bestTune$lambda # 5

rr_14k  <- linearRidge(off_survive ~ ., data = all_data_RRS, lambda = lambda_14k)
rr_all  <- linearRidge(off_survive ~ ., data = combine_data_RRS, lambda = lambda_all)
rr_meth <- linearRidge(off_survive ~ ., data = methyl_only, lambda = lambda_meth)
rr_gen  <- linearRidge(off_survive ~ ., data = gen_only, lambda = lambda_gen)
rr_env  <- linearRidge(off_survive ~ ., data = env_only, lambda = lambda_env)

R2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

y <- combine_data_RRS$off_survive

pred_14k  <- predict(rr_14k)
pred_all  <- predict(rr_all)
pred_meth <- predict(rr_meth)
pred_gen  <- predict(rr_gen)
pred_env  <- predict(rr_env)

r2_14k  <- R2(y, pred_14k) 
r2_all  <- R2(y, pred_all)
r2_meth <- R2(y, pred_meth)
r2_gen  <- R2(y, pred_gen)
r2_env  <- R2(y, pred_env)

r2_14k; r2_all; r2_meth; r2_gen; r2_env

data.frame(
  Model = c("All", "Methylation", "Genetic", "Environment"),
  R2 = c(r2_all, r2_meth, r2_gen, r2_env)
)

r2_14k # 0.6586617

# Model          R2
# 1         All 0.167349865
# 2 Methylation 0.145559521
# 3     Genetic 0.003832832
# 4 Environment 0.136794832

# CLEAR
remove(methyl_only, gen_only, env_only, combine_data_RRS, ridge_train, lambda_all, lambda_env, lambda_gen, lambda_meth, 
       rr_all, rr_meth, rr_env, rr_gen, r2_all, r2_meth, r2_gen, r2_env)

## CODE FOR LONGEVITY
# Data subsets
methyl_only <- data.frame(methyl_subset, longevity = sig_methyl_data$longevity)
gen_only <- data.frame(sig_methyl_data[, c("lineage", "pop", "longevity")])
env_only <- data.frame(sig_methyl_data[, c("month", "year", "rain", "temp", "totpopsz", "longevity")])

# Combined (already defined)
combine_data_long <- data.frame(methyl_subset, sig_methyl_data[,c("lineage", "pop", "totpopsz",
                                                                 "month", "year", "rain", "temp",
                                                                 "longevity")])

all_data_long <- data.frame(winds, meta[,c("lineage", "pop", "totpopsz",
                                          "month", "year", "rain", "temp",
                                          "longevity")])

library(caret)
library(ridge)

ridge_train <- function(df) {
  caret::train(longevity ~ .,
               data = df,
               method = "glmnet",
               tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(2, 3, length = 100)),
               trControl = trainControl(method = "repeatedcv", number = 10, repeats = 100))
}

long_all  <- ridge_train(combine_data_long)
long_meth <- ridge_train(methyl_only)

ridge_train <- function(df) {
  caret::train(longevity ~ .,
               data = df,
               method = "glmnet",
               tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(0, 1, length = 100)),
               trControl = trainControl(method = "repeatedcv", number = 10, repeats = 100))
}

long_env  <- ridge_train(env_only)
long_14k  <- ridge_train(all_data_long)

ridge_train <- function(df) {
  caret::train(longevity ~ .,
               data = df,
               method = "glmnet",
               tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(1, 2, length = 100)),
               trControl = trainControl(method = "repeatedcv", number = 10, repeats = 100))
}
long_gen  <- ridge_train(gen_only)

lambda_14k  <- long_14k$bestTune$lambda # 6
lambda_all  <- long_all$bestTune$lambda # 129
lambda_meth <- long_meth$bestTune$lambda # 120 -- may be an error here -- double check
lambda_gen  <- long_gen$bestTune$lambda # 91 -- may be an effor here -- double check
lambda_env  <- long_env$bestTune$lambda # 2

rr_14k  <- linearRidge(longevity ~ ., data = all_data_long, lambda = lambda_14k)
rr_all  <- linearRidge(longevity ~ ., data = combine_data_long, lambda = lambda_all)
rr_meth <- linearRidge(longevity ~ ., data = methyl_only, lambda = lambda_meth)
rr_gen  <- linearRidge(longevity ~ ., data = gen_only, lambda = lambda_gen)
rr_env  <- linearRidge(longevity ~ ., data = env_only, lambda = lambda_env)

R2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

y <- combine_data_long$longevity

pred_14k  <- predict(rr_14k)
pred_all  <- predict(rr_all)
pred_meth <- predict(rr_meth)
pred_gen  <- predict(rr_gen)
pred_env  <- predict(rr_env)

r2_14k  <- R2(y, pred_14k)
r2_all  <- R2(y, pred_all)
r2_meth <- R2(y, pred_meth)
r2_gen  <- R2(y, pred_gen)
r2_env  <- R2(y, pred_env)

r2_all; r2_meth; r2_gen; r2_env

unique_env  <- r2_env
unique_gen  <- r2_gen - r2_env
unique_meth <- r2_all - r2_gen
shared_total <- r2_all - (unique_env + unique_gen + unique_meth)

data.frame(
  Model = c("All", "Methylation", "Genetic", "Environment"),
  R2 = c(r2_all, r2_meth, r2_gen, r2_env)
)

r2_14k # 0.9985887

# Model           R2
# 1         All 0.2539449135
# 2 Methylation 0.2617248998
# 3     Genetic 0.0008576057
# 4 Environment 0.0894603485