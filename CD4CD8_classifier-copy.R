#ML model is loaded here, applied to car t data in donor_analysis.R

library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)
library(Seurat)
library(xgboost)
library(SHAPforxgboost)
library(caret)

#Plot properly
options(bitmapType="cairo") 

# Load Seurat object for training model
tcell <- LoadSeuratRds("/pfs/stor10/users/home/m/mb223/mystore/cartdata/tcelldata.ccreg.Rdata")

# Visualize gene markers
FeaturePlot(tcell, features = c("CD4", "CD8A"))

# Extract expression data for CD4 and CD8A genes
df <- FetchData(tcell, vars=c("CD4","CD8A"))

# Assign known CD4/CD8 cell types
df$tcell_type <- NA
df$tcell_type[df$CD4>0 & df$CD8A==0] <- 0
df$tcell_type[df$CD8A>0 & df$CD4==0] <- 1

# Fetch gene expression data, excluding CD4/CD8
toc <- FetchData(tcell, vars = setdiff(VariableFeatures(tcell), c("CD4", "CD8")))
toc <- as.matrix(toc)

# Create training and test sets from known CD4/CD8 cells
known_idx <- which(!is.na(df$tcell_type))
tosplit <- round(runif(length(known_idx)))
train_idx <- known_idx[tosplit==0]
test_idx <- known_idx[tosplit==1]

# Prepare XGBoost input matrices
X_train <- toc[train_idx,]
y_train <- df$tcell_type[train_idx]
X_test <- toc[test_idx,]

train_input <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
test_input <- xgb.DMatrix(data = as.matrix(X_test))

# Set XGBoost parameters
params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  eta = 0.3,
  max_depth = 4
)

# Train XGBoost model
CDClassifier <- xgboost(data = train_input, 
                     params = params, 
                     nrounds = 100, 
                     verbose = 1)

# Predict on test data
predictions <- predict(CDClassifier, test_input)
preds_outcome <- ifelse(predictions < 0.5, 0, 1)

df$tcell_pred <- NA
df$tcell_pred[test_idx] <- preds_outcome

# Visualize predictions 
table(df$tcell_type[test_idx], df$tcell_pred[test_idx], useNA = "ifany")

# Compute confusion matrix and accuracy 
conf_mat <- table(df$tcell_type[test_idx], df$tcell_pred[test_idx], useNA = "ifany")
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
print(accuracy)  

# Predict labels for all cells
X_all <- toc
all_input <- xgb.DMatrix(data = as.matrix(X_all))
preds <- predict(CDClassifier, all_input)
final_preds <- ifelse(preds < 0.5, 0, 1)
df$tcell_pred <- final_preds
table(df$tcell_type, df$tcell_pred, useNA = "ifany")


# Classification model evaluation 
xgbcv <- xgb.cv( params = params, 
                 data = train_input, 
                 nrounds = 100, 
                 nfold = 5, 
                 showsd = TRUE, 
                 stratified = TRUE, 
                 print_every_n = 10, 
                 early_stop_round = 20, 
                 maximize = FALSE)

best_nrounds <- which.min(xgbcv$evaluation_log$test_logloss_mean)
best_mlogloss <- min(xgbcv$evaluation_log$test_logloss_mean)
print(c(best_nrounds, best_mlogloss))

# Plot cross-validation results
cv_results <- xgbcv$evaluation_log
cv_long <- reshape2::melt(cv_results, id.vars = "iter")

ggplot(cv_long, aes(x = iter, y = value, color = variable)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = best_nrounds, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "XGBoost Cross-Validation Log-Loss",
       x = "Iteration",
       y = "Log-Loss",
       color = "Metric") +
  scale_color_manual(values = c("blue", "lightblue", "red", "pink"))


# SHAP analysis for feature importance
shap_values <- shap.values(CDClassifier, X_train)

pdf(file = "./mystore/cartdata/Plots/CDShap.pdf",  width = 6,  height = 5)
shap.plot.summary.wrap1(CDClassifier, X_train, top_n = 20) +
  ggtitle("SHAP Summary Plot for CD4/CD8 Classifier")
dev.off()
# ggsave("figs/CD4_CD8_summary.png", plot = last_plot(), width = 6, height = 4)


# Add classification results to Seurat metadata
tcell.ccreg <- AddMetaData(tcell.ccreg, metadata = df$tcell_pred, col.name = "CD4CD8_Pred")

# Save updated object
SaveSeuratRds(tcell.ccreg, file = "/pfs/stor10/users/home/a/awallin/ondemand/wallin/tcell_pred1.Rdata")

#


