#ML model is loaded here, applied to car t data in donor_analysis.R

library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)
library(Seurat)
library(xgboost)
library(SHAPforxgboost)

#Plot properly
options(bitmapType="cairo") 

# Load Seurat object containing predicted T cell types
tcell <- LoadSeuratRds("/pfs/stor10/users/home/m/mb223/mystore/cartdata/tcell_pred1.Rdata")

# Fetch CAR data
df2 <- FetchData(tcell, vars = "detect_CAR") 
df2$detect_CAR[df2$detect_CAR==FALSE] <- 0
df2$detect_CAR[df2$detect_CAR==TRUE] <- 1

# Fetch gene expression data excluding CAR genes
all_features <- VariableFeatures(tcell)
CAR_genes <- c("virus-M28z", "virus-MBBz", "EGFR")  
filtered_features <- setdiff(all_features, CAR_genes)

# Extract expression matrix
emat <- FetchData(tcell, vars = filtered_features)
emat <- as.matrix(emat)

# Identify known CAR labels and split data into train/test
car_index <- which(!is.na(df2$detect_CAR))
set.seed(42)
train_index <- sample(car_index, size = round(0.8 * length(car_index)))
test_index <- setdiff(car_index, train_index)

# Prepare train/test sets
X_training <- emat[train_index,]
y_training <- df2$detect_CAR[train_index]
X_testing <- emat[test_index,]
y_testing <- df2$detect_CAR[test_index]

train_mat <- xgb.DMatrix(data = as.matrix(X_training), label = y_training)
test_mat <- xgb.DMatrix(data = as.matrix(X_testing), label = y_testing)

# Define XGBoost parameters
parameters <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  eta = 0.3,  
  max_depth = 4
)

# Train XGBoost model
car_classifier <- xgboost(data = train_mat, 
                          params = parameters, 
                          nrounds = 100, 
                          verbose = 1)

# Preciction evaluation
pred <- predict(car_classifier, test_mat)
pred_outcome <- ifelse(pred < 0.5, 0, 1)
df2$CAR_pred <- NA
df2$CAR_pred[test_index] <- pred_outcome

# Compute confustion matrix and accuracy 
conf_mat <- table(df2$detect_CAR[test_index], df2$CAR_pred[test_index], useNA = "ifany")
print(conf_mat)
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
print(paste("Model Accuracy:", round(accuracy, 4)))

# Predict CAR status for all cells
all_input <- xgb.DMatrix(data = as.matrix(emat))
all_preds <- predict(car_classifier, all_input)
final_preds <- ifelse(all_preds < 0.5, 0, 1)
df2$CAR_pred <- final_preds

# Summarize classification results
table(df2$detect_CAR, df2$CAR_pred, useNA = "ifany")

# Model evaluation
cv <- xgb.cv(params = parameters, 
             data = train_mat, 
             nrounds = 100, 
             nfold = 5, 
             showsd = TRUE, 
             stratified = TRUE, 
             print_every_n = 10, 
             early_stop_round = 20, 
             maximize = FALSE)

n_rounds <- which.min(cv$evaluation_log$test_logloss_mean)
eval_metrics <- min(cv$evaluation_log$test_logloss_mean)
print(c(n_rounds, eval_metrics))

# Plot cross-validation performance
cv_eval <- cv$evaluation_log
cv_formatted <- reshape2::melt(cv_eval, id.vars = "iter")

ggplot(cv_formatted, aes(x = iter, y = value, color = variable)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = n_rounds, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "XGBoost Cross-Validation Log-Loss",
       x = "Iteration",
       y = "Log-Loss",
       color = "Metric") +
  scale_color_manual(values = c("blue", "lightblue", "red", "pink"))

# SHAP analysis for feature importance
shap_values <- shap.values(car_classifier, X_training)
shap.plot.summary.wrap1(car_classifier, X_training, top_n = 20) +
  ggtitle("SHAP Summary Plot for Car Classifier")
# ggsave("figs/CAR_summary.png", plot = last_plot(), width = 5, height = 4)

# Add CAR classification results to Seurat metadata
tcell <- AddMetaData(tcell, metadata = df2$CAR_pred, col.name = "CAR_Pred")
# head(tcell@meta.data)

# Save updated Seurat object
SaveSeuratRds(tcell, file = "/pfs/stor10/users/home/a/awallin/ondemand/wallin/tcell_pred2.Rdata")
