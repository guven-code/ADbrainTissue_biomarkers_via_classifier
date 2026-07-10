# =====================================================================
# 1. SETUP & LIBRARIES
# =====================================================================
message("Loading Libraries...")
required_pkgs <- c(
  "GEOquery", "caret", "glmnet", "org.Hs.eg.db", "AnnotationDbi",
  "pROC", "data.table", "edgeR", "sva", "PRROC", "ggplot2",
  "limma", "reshape2", "pheatmap", "gridExtra", "dplyr",
  "viridis", "tidyr", "clusterProfiler", "stringr", "ggrepel",
  "ComplexHeatmap", "circlize"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GEOquery","org.Hs.eg.db","AnnotationDbi","clusterProfiler","edgeR","sva")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE)
    } else install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
library(dplyr)
library(tidyr)
library(ggplot2)
library(caret)
library(purrr)


# =====================================================================
# 2. UTILITIES
# =====================================================================
get_counts <- function(gse_id, file_name) {
  base <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
  path <- paste(base, paste0("acc=", gse_id), paste0("file=", file_name), sep="&")
  as.matrix(data.table::fread(path, header=TRUE), rownames=1)
}

get_label <- function(ph, pattern_col, pos_pattern) {
  col <- grep(pattern_col, colnames(ph), ignore.case=TRUE, value=TRUE)[1]
  raw <- as.character(ph[[col]])
  factor(ifelse(grepl(pos_pattern, raw, ignore.case=TRUE), "AD", "Control"),
         levels=c("Control","AD"))
}

filter_counts <- function(x) {
  x[rowSums(x >= 10) >= max(2, round(ncol(x)*0.1)), ]
}

prep_matrix <- function(mat, genes) {
  m <- mat[genes, , drop=FALSE]
  t(m)
}

# =====================================================================
# 3. LOAD DATA
# =====================================================================
gse1 <- getGEO("GSE125583", GSEMatrix=TRUE)[[1]]
ph1  <- pData(gse1)
tbl1 <- get_counts("GSE125583","GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz")

gse2 <- getGEO("GSE153873", GSEMatrix=TRUE)[[1]]
ph2  <- pData(gse2)
tbl2 <- get_counts("GSE153873","GSE153873_raw_counts_GRCh38.p13_NCBI.tsv.gz")

ph1 <- ph1[-38,]
tbl1 <- tbl1[,-38]

y1 <- get_label(ph1, "title|disease", "Alzheimer|AD")
y2 <- factor(ifelse(grepl("Alzheimer",
                          ph2[[grep("disease state", colnames(ph2), value=TRUE)[1]]]),
                    "AD","Control"),
             levels=c("Control","AD"))

age1 <- as.numeric(ph1$`age:ch1`)
sex1 <- factor(ifelse(ph1$`Sex:ch1` %in% c("F","Female"), "Female","Male"))

# =====================================================================
# 4. NORMALIZATION (NO TEST INFO USED IN TRAINING PIPELINE)
# =====================================================================
set1 <- edgeR::cpm(filter_counts(tbl1), log=TRUE, prior.count=1)
set2 <- edgeR::cpm(filter_counts(tbl2), log=TRUE, prior.count=1)
set1[1:3,1:3]; set2[1:3,1:3];
# IMPORTANT FIX:
# Use ONLY training dataset to define feature space
#genes_all <- rownames(set1)

# Map symbols safely (training-only definition)
symbols1 <- mapIds(org.Hs.eg.db,
                   keys=rownames(set1),
                   column="SYMBOL",
                   keytype="ENTREZID",
                   multiVals="first")

symbols2 <- mapIds(org.Hs.eg.db,
                   keys=rownames(set2),
                   column="SYMBOL",
                   keytype="ENTREZID",
                   multiVals="first")

remove_genes1 <- grepl("^MT-|^RPL|^RPS|^SNORD|^SNORA|^RNU|^RNA", symbols1)
remove_genes2 <- grepl("^MT-|^RPL|^RPS|^SNORD|^SNORA|^RNU|^RNA", symbols2)

set1 <- set1[!remove_genes1, ]
symbols1 <- symbols1[!remove_genes1]

set2 <- set2[!remove_genes2, ]
symbols2 <- symbols2[!remove_genes2]


  
keep1 <- !is.na(symbols1)
set1 <- set1[keep1, ]
symbols1 <- symbols1[keep1]

keep2 <- !is.na(symbols2)
set2 <- set2[keep2, ]
symbols2 <- symbols2[keep2]

  
set1 <- set1[!duplicated(symbols1), ]
symbols1 <- symbols1[!duplicated(symbols1)]

set2 <- set2[!duplicated(symbols2), ]
symbols2 <- symbols2[!duplicated(symbols2)]

  
rownames(set1) <- symbols1
rownames(set2) <- symbols2

common_genes <- intersect(rownames(set1), rownames(set2))

set1 <- set1[common_genes, ]
set2 <- set2[common_genes, ]


  
sum(is.na(rownames(set1)))      # should be 0
sum(duplicated(rownames(set1))) # should be 0

sum(is.na(rownames(set2)))      # should be 0
sum(duplicated(rownames(set2))) # should be 0



X1 <- t(set1)
X2 <- t(set2)
dim(X2)

# =====================================================================
# 5. TRAIN/TEST SPLIT (NO FEATURE SELECTION BEFORE THIS)
# =====================================================================
set.seed(1234)
idx <- createDataPartition(y1, p=0.8, list=FALSE)

train_x_raw <- X1[idx,]
train_y     <- y1[idx]
train_age   <- age1[idx]
train_sex   <- sex1[idx]

test_x_raw  <- X1[-idx,]
test_y      <- y1[-idx]

# External dataset untouched
ext_x_raw <- X2
ext_y     <- y2

train_x_raw[1:3,1:3]
summary(rowMeans(t(train_x_raw)))



# =====================================================================
# 6. NORMALIZATION (TRAIN ONLY FIT)
# =====================================================================
library(caret)
pp <- preProcess(train_x_raw, method=c("center","scale"))

train_x <- predict(pp, train_x_raw)
test_x  <- predict(pp, test_x_raw)
ext_x   <- predict(pp, ext_x_raw)

# =====================================================================
# 7. FEATURE SELECTION (IMPORTANT FIX: INSIDE CV ONLY)
# =====================================================================

y <- factor(train_y, levels=c("Control","AD"))

design <- model.matrix(~ y + train_age + train_sex)

fit <- limma::lmFit(t(train_x), design)
fit <- limma::eBayes(fit)

tt <- limma::topTable(fit, coef="yAD", number=Inf)

sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.8, ]

up   <- rownames(sig)[sig$logFC > 0.8]
down <- rownames(sig)[sig$logFC < -0.8]

selected_genes<-c(head(up,15), head(down,15))
# IMPORTANT: we will NOT precompute genes globally anymore
# feature selection will be embedded in modeling stage conceptually

# =====================================================================
# 8. MODELING WITH FIXED FEATURE SPACE APPROACH (SAFE VERSION)
# =====================================================================

train_df <- data.frame(train_x[, selected_genes], condition = train_y)
test_df  <- data.frame(test_x[, selected_genes], condition = test_y)
ext_df   <- data.frame(ext_x[, selected_genes], condition = ext_y)


# =====================================================================
# 6. MODEL TRAINING (Fixes Critique 1: Removed sampling="down" leak)
# =====================================================================
# sampling=NULL completely safeguards fold isolation. Class weights adjust for imbalances.
ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, 
                     summaryFunction = twoClassSummary, sampling = NULL, savePredictions = "final") 

# Generate explicit model-specific observation weights based on the training split distribution
weights_vec <- ifelse(train_df$condition == "AD", 1 / table(train_df$condition)["AD"], 1 / table(train_df$condition)["Control"])

grid_glmnet <- expand.grid(alpha = seq(0, 1, by = 0.2), lambda = seq(0.001, 0.1, length = 10))
grid_gbm    <- expand.grid(n.trees = c(50, 100, 150), interaction.depth = c(1, 3, 5), shrinkage = c(0.01, 0.1), n.minobsinnode = 10)
grid_nnet   <- expand.grid(size = c(1, 3, 5, 10), decay = c(0.1, 0.5, 1.0))

# Note: Weights are handled natively by GLMNET, GBM, and NNET; caret silently drops them for SVM, PLS, and RF.
set.seed(1234); fit_glmnet <- train(condition ~ ., train_df, method="glmnet", trControl=ctrl, metric="ROC", tuneGrid = grid_glmnet, weights = weights_vec)
set.seed(1234); fit_svm    <- train(condition ~ ., train_df, method="svmRadial", trControl=ctrl, metric="ROC", tuneLength = 10, weights = weights_vec) 
set.seed(1234); fit_pls    <- train(condition ~ ., train_df, method="pls", trControl=ctrl, metric="ROC", tuneLength = 15, weights = weights_vec)
#set.seed(1234); fit_rf     <- train(condition ~ ., train_df, method="rf", trControl=ctrl, metric="ROC", tuneLength = 5, weights = weights_vec)
set.seed(1234); fit_gbm    <- train(condition ~ ., train_df, method="gbm", trControl=ctrl, metric="ROC", tuneGrid = grid_gbm, verbose=FALSE, weights = weights_vec)
set.seed(1234); fit_nnet   <- train(condition ~ ., data = train_df, method = "nnet", trControl = ctrl, metric = "ROC", tuneGrid = grid_nnet, trace = FALSE, maxit = 500, weights = weights_vec)

models <- list(GLMNET = fit_glmnet, SVM = fit_svm, PLS = fit_pls,
               #RF = fit_rf,
               GBM = fit_gbm, NNET = fit_nnet)
datasets <- list(TRAIN=train_df, TEST=test_df, EXT=ext_df) # EXT2=set3_df, EXT3=set4_df)

# =====================================================================
# 7. METRICS SUMMARY AND PURIFIED THRESHOLDS
# =====================================================================
res_samples <- caret::resamples(models)
res_summary <- summary(res_samples)

res_samples <- caret::resamples(models)

# Isolate the exact columns containing the ROC metric from the raw values matrix
roc_cols <- res_samples$values[, grep("~ROC$", colnames(res_samples$values)), drop = FALSE]

# Strip out the "~ROC" suffix from column names to get clean model names (e.g., "GLMNET")
clean_names <- gsub("~ROC$", "", colnames(roc_cols))

# Compute Means and SDs directly from the raw fold data to bypass summary table limits
cv_means <- colMeans(roc_cols, na.rm = TRUE)
cv_sds   <- apply(roc_cols, 2, sd, na.rm = TRUE)

# Reassign names to match your models list exactly
names(cv_means) <- clean_names
names(cv_sds)   <- clean_names
# Fix: Generate un-contaminated training thresholds strictly via CV fold predictions (Eliminates Optimistic Training Bias)
frozen_thresholds <- lapply(models, function(m) {
  pred_df <- m$pred
  best_params <- m$bestTune
  
  # Strictly isolate rows matching the optimized hyperparameter set
  for(param in names(best_params)){
    pred_df <- pred_df[pred_df[[param]] == best_params[[param]], ]
  }
  
  # Calculate threshold based on cross-validated holdout folds
  roc_obj <- pROC::roc(pred_df$obs, pred_df$AD, quiet=TRUE, levels=c("Control", "AD"), direction="<")
  as.numeric(pROC::coords(roc_obj, "best", ret="threshold", best.method="youden")[1])
})

table_data <- list()
for (d in names(datasets)) {
  for (m in names(models)) {
    data_curr <- datasets[[d]]
    truth     <- factor(data_curr$condition, levels = c("Control", "AD"))
    probs     <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    
    preds     <- factor(ifelse(probs >= frozen_thresholds[[m]], "AD", "Control"), levels = c("Control", "AD"))
    cm        <- confusionMatrix(preds, truth, positive = "AD")
    
    table_data[[paste(d, m, sep = "_")]] <- data.frame(
      Model = m, Dataset = d, Accuracy = round(cm$overall["Accuracy"] * 100, 1),
      Precision = round(cm$byClass["Pos Pred Value"] * 100, 1), Recall = round(cm$byClass["Sensitivity"] * 100, 1),
      F1_score = round(cm$byClass["F1"] * 100, 1), Mean_CV = round(cv_means[m] * 100, 1), SD_CV = round(cv_sds[m], 3)
    )
  }
}
final_output_table <- do.call(rbind, table_data)

# Convert the frozen_thresholds list into a clean data frame
thresholds_df <- data.frame(
  Model = names(frozen_thresholds),
  Optimal_Threshold = unlist(frozen_thresholds),
  row.names = NULL
)

# Print the table to the console
print("--- Optimized Thresholds (Youden Index from CV) ---")
print(thresholds_df)


library(pROC)

auc_data <- list()

for (d in names(datasets)) {
  for (m in names(models)) {
    data_curr <- datasets[[d]]
    truth      <- factor(data_curr$condition, levels = c("Control", "AD"))
    
    # Get predicted probabilities for the "AD" class
    probs      <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    
    # Calculate the AUC using pROC
    roc_obj    <- pROC::roc(truth, probs, quiet = TRUE, levels = c("Control", "AD"), direction = "<")
    calculated_auc <- as.numeric(pROC::auc(roc_obj))
    
    # If it's the TRAIN dataset, we can also grab the Cross-Validation AUC for comparison
    cv_auc_val <- ifelse(d == "TRAIN", cv_means[m], NA)
    
    auc_data[[paste(d, m, sep = "_")]] <- data.frame(
      Model = m,
      Dataset = d,
      Threshold = round(frozen_thresholds[[m]], 3),
      Dataset_AUC = round(calculated_auc, 3),
      CV_Holdout_AUC = if(!is.na(cv_auc_val)) round(cv_auc_val, 3) else "-"
    )
  }
}

# Combine and print the results
auc_summary_table <- do.call(rbind, auc_data)
rownames(auc_summary_table) <- NULL

print("--- AUC and Threshold Summary ---")
print(auc_summary_table)

library(pROC)
library(caret)

table_data <- list()

for (d in names(datasets)) {
  for (m in names(models)) {
    data_curr <- datasets[[d]]
    truth     <- factor(data_curr$condition, levels = c("Control", "AD"))
    
    # 1. Get predicted probabilities and class assignments
    probs     <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    preds     <- factor(ifelse(probs >= frozen_thresholds[[m]], "AD", "Control"), levels = c("Control", "AD"))
    
    # 2. Confusion Matrix & Accuracy Confidence Intervals
    cm        <- confusionMatrix(preds, truth, positive = "AD")
    acc       <- round(cm$overall["Accuracy"], 2)
    acc_lower <- round(cm$overall["AccuracyLower"], 2)
    acc_upper <- round(cm$overall["AccuracyUpper"], 2)
    
    # Format Accuracy with its CI
    acc_with_ci <- paste0(acc, " (95% CI ", acc_lower, "–", acc_upper, ")")
    
    # Other metrics
    precision <- round(cm$byClass["Pos Pred Value"] * 100, 1)
    recall    <- round(cm$byClass["Sensitivity"] * 100, 1)
    f1_score  <- round(cm$byClass["F1"] * 100, 1)
    
    # 3. Calculate AUC and AUC Confidence Intervals
    roc_obj   <- pROC::roc(truth, probs, quiet = TRUE, levels = c("Control", "AD"), direction = "<")
    auc_val   <- as.numeric(pROC::auc(roc_obj))
    
    # pROC returns a vector of length 3: [Lower bound, Median/AUC, Upper bound]
    auc_ci    <- pROC::ci.auc(roc_obj)
    
    # Format AUC with its CI
    auc_with_ci <- paste0(round(auc_val, 2), " (95% CI ", round(auc_ci[1], 2), "–", round(auc_ci[3], 2), ")")
    
    # 4. Store all metrics cleanly
    table_data[[paste(d, m, sep = "_")]] <- data.frame(
      Model = m, 
      Dataset = d, 
      Accuracy_CI = acc_with_ci,
      Precision = precision, 
      Recall = recall, 
      F1_score = f1_score, 
      AUC_CI = auc_with_ci,
      Threshold = round(frozen_thresholds[[m]], 3),
      row.names = NULL
    )
  }
}

final_output_table_with_ci <- do.call(rbind, table_data)
rownames(final_output_table_with_ci) <- NULL

print("--- Final Model Performance Table with 95% Confidence Intervals ---")
print(final_output_table_with_ci[, c("Model", "Dataset", "Accuracy_CI", "AUC_CI", "Threshold")])


library(MLeval)

# A) Evaluate your individual models on the TEST set
# This creates ROC, PR, and Calibration curves automatically
tiff("reanalyze_roc_ModelsSET1.tiff",width=6,height=6, unit= "in", res = 300 )
res <- evalm(list(fit_glmnet,fit_svm,fit_pls,fit_gbm,fit_nnet), 
             gnames = c('GLMNET', 'SVM','PLS','GBM','NNET'),
             title = "Multi-Model Evaluation: Training Cohort")
dev.off()


models_compare <- caret::resamples(list(
  GLMNET = fit_glmnet,
  SVM    = fit_svm,
  PLS    = fit_pls,
  #RF     = fit_rf,
  GBM    = fit_gbm,
  NNET   = fit_nnet
))
print(summary(models_compare))

# 2. Define High-Legibility Graphical Settings
# Adjusting lattice settings for publication standards
par.settings <- list(
  axis.text = list(cex = 1.4, font = 2),       # Bold and Large Y-axis (Model Names)
  axis.title = list(cex = 1.5, font = 2),      # Bold and Large Axis Titles
  strip.text = list(cex = 1.4, font = 2),      # Bold and Large Panel Titles (ROC, Sens, Spec)
  plot.symbol = list(col = "black", pch = 16), # Clean black dots for outliers
  box.rectangle = list(fill = "steelblue", alpha = 0.5) # Professional blue boxes
)

# 3. Define Plotting Scales
# We keep relation="free" so each metric (ROC, Sens, Spec) has its own range
scales <- list(
  x = list(relation = "free", cex = 1.2, font = 2), 
  y = list(relation = "free", cex = 1.2, font = 2)
)
# Plot with the modified settings (Removed the conflicting xlab)
lattice::bwplot(models_compare, 
                scales = scales, 
                par.settings = par.settings,
                main = list(label = "Cross-Validation Performance Summary", cex = 1.6, font = 2))
tiff(file="reanalyze_Figure_Model_Comparison_Yaxis.tiff", 
     unit="in", res=300, width=10, height=7, compression="lzw")

# 1. Create the plot object without the xlab error
p <- bwplot(models_compare, scales = scales, par.settings = par.settings)

# 2. Update the labels and title manually
update(p, 
       main = list(cex = 1.6, font = 2, label = "Cross-Validation Performance Summary"),
       xlab = list(cex = 1.5, font = 2, label = "Performance Metric Value"))
dev.off()

#####DEGS volcano plot:
library(ggplot2)
library(ggrepel)

# 1. Create a category column in your main results table 'tt'
# --- A. VOLCANO PLOT (Upgraded to High-Contrast Blueprint) ---
y <- factor(train_y, levels=c("Control","AD"))

design <- model.matrix(~ y + train_age + train_sex)

fit <- limma::lmFit(t(train_x), design)
fit <- limma::eBayes(fit)

tt <- limma::topTable(fit, coef="yAD", number=Inf)

sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.8, ]

up   <- rownames(sig)[sig$logFC > 0.8]
down <- rownames(sig)[sig$logFC < -0.8]

selected_genes_30<-c(head(up,15), head(down,15))
# 1. Update category labeling for optimal contrast mapping
tt$category <- "Not Significant"
tt$category[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.8] <- "Significant (Other)"
tt$category[rownames(tt) %in% selected_genes_30] <- "Selected 30 Biomarkers"

# 2. Define specific visual weights to force background data to recede
tt$point_size <- 0.5
tt$point_size[tt$category == "Significant (Other)"] <- 0.8
tt$point_size[tt$category == "Selected 30 Biomarkers"] <- 2.5

# Uncomment the following lines if you want to export the figure directly to TIFF
tiff("reanalyze_Volcano_HighContrast.tiff", width=8, height=6, units="in", res=300)

# 3. Generate high-contrast plot layout
volcano_p <- ggplot(tt, aes(x=logFC, y=-log10(adj.P.Val), color=category, size=point_size)) +
  # Apply the variable point sizes defined dynamically in the dataframe
  geom_point(alpha=0.6) + 
  scale_size_identity() +
  
  # Ensure the 30 biomarkers are labeled clearly without clipping or crowded overlaps
  geom_text_repel(data=subset(tt, category=="Selected 30 Biomarkers"), 
                  aes(label=rownames(subset(tt, category=="Selected 30 Biomarkers"))),
                  size=3.2, 
                  color="black", 
                  fontface="bold",
                  box.padding = 0.6, 
                  max.overlaps = Inf,
                  family="sans") +
  
  # Inject higher contrast color mapping palette
  scale_color_manual(values=c("Selected 30 Biomarkers"="red3", 
                              "Significant (Other)"="deepskyblue3", 
                              "Not Significant"="grey85")) +
  
  # Zoom in on the X-axis window to separate vertical gene clusters beautifully
  coord_cartesian(xlim = c(-1.5, 1.5)) + 
  
  theme_bw() + 
  theme(
    text = element_text(family="sans"),
    plot.title = element_text(hjust = 0.5, face="bold", size=14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(title="Volcano Plot: Distribution of 30 Selected Biomarkers",
       x="log2 Fold Change", 
       y="-log10 Adjusted P-value",
       color="Gene Category") +
  
  # Add signature statistical cutoff reference lines
  geom_vline(xintercept=c(-0.8, 0.8), linetype="dotted", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color="grey50")

# Render the layout canvas
print(volcano_p)

dev.off() # Uncom


###32.nd option volcano:
# 1. Update the categories to distinguish direction
tt$category <- "Not Significant"
tt$category[tt$adj.P.Val < 0.05 & tt$logFC > 0.8] <- "Up-regulated"
tt$category[tt$adj.P.Val < 0.05 & tt$logFC < -0.8] <- "Down-regulated"

# -----------------------------
# 1. Prepare data
# -----------------------------
volcano_df <- tt %>%
  mutate(
    gene = rownames(tt),
    negLogP = -log10(adj.P.Val),
    status = case_when(
      adj.P.Val < 0.05 & logFC > 0.8  ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.8 ~ "Down",
      TRUE ~ "NotSig"
    ),
    highlight30 = ifelse(gene %in% selected_genes_30, "Top30", "Other")
  )
#sum(tt$adj.P.Val <0.05 & tt$logFC > 0.8 ) 209 up-genes
#sum(tt$adj.P.Val <0.05 & tt$logFC  < -0.8 ) 436 down-genes
# -----------------------------
# 2. Volcano Plot
# -----------------------------
tiff("reanalyze_Volcano_HighContrast.tiff", width=8, height=6, units="in", res=300)
library(ggplot2)
library(ggrepel)
library(dplyr)

ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  
  # Base points (all genes) - Smaller size (0.8)
  geom_point(aes(color = status), alpha = 0.6, size = 0.8) +
  
  # Highlight top 30 genes - Mapped inside aes() and smaller size (1.8)
  geom_point(
    data = volcano_df %>% filter(highlight30 == "Top30"),
    aes(color = "Biomarker genes"),   # <-- Moved inside aes() to fix the legend
    size = 1.8
  ) +
  
  # Labels for top 30 genes
  geom_text_repel(
    data = volcano_df %>% filter(highlight30 == "Top30"),
    aes(label = gene),
    size = 3,                         # <-- Slightly smaller text to match smaller dots
    box.padding = 0.4,
    max.overlaps = 100
  ) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # Colors (Now "Biomarker genes" will map correctly)
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NotSig" = "grey",
    "Biomarker genes" = "purple"
  )) +
  
  # Theme adjustments
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  
  # Labels
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Expression",
    subtitle = paste0("Number of up-regulated genes is ", 209, 
                      "\nNumber of down-regulated genes is ", 436)
  )
dev.off()


####re-analyze heatmap

library(ComplexHeatmap)
library(circlize)

# 1. Prepare the data matrix (30 selected genes)
# Assuming train_x and train_y are sorted or ordered as in your pheatmap code
ordered_indices <- order(train_y)
plot_mat <- t(train_x[ordered_indices, selected_genes_30])

# 2. Define the Color Map for the expression values (Standard Z-score range)
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# 3. Define the Top Annotation (The "Condition" bar)
# Matching your requested yellow/black scheme
top_ann = HeatmapAnnotation(
  Condition = train_y[ordered_indices],
  col = list(Condition = c("Control" = "yellow", "AD" = "black")),
  show_legend = TRUE
)

tiff("reanalyze_Figure3_Heatmap_Complex.tiff", 
     width = 10,
     height = 8,
     units = "in",
     res = 300,
     compression = "lzw")
# # 4. Generate the Heatmap
ch<-ComplexHeatmap::Heatmap(
  plot_mat, 
  name = "Z-score", 
  top_annotation = top_ann,
  col = col_fun,
  
  # Clustering & Display
  cluster_columns = FALSE,          # Keep your custom ordering
  cluster_rows = TRUE,             # Group similar genes together
  show_column_names = FALSE,       # Often too crowded for publications
  show_row_names = TRUE,
  
  # Aesthetics
  row_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
  column_title = "Differential Expression Profile: Top Predictive Genes",
  column_title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial"),
  
  # Legend parameters
  heatmap_legend_param = list(title = "Rel. Expression",fontsize=12)
)
draw(ch)
dev.off()

####GBM impotance
# =====================================================================
# 8. BIOLOGICAL INTERPRETATION (GBM & ENRICHMENT)
# =====================================================================
# GBM Variable Importance
library(gbm)
gbm_imp <- varImp(models$GBM, scale=TRUE)
tiff("GBM_Importance.tiff", width=8, height=6, units="in", res=300)
gbm_imp <- varImp(models$GBM, scale=TRUE)
# Store the plot first
p <- plot(gbm_imp, top = 20, main = "Top 20 Important Genes (GBM)")

# Update with larger text sizes
update(p, 
       xlab = list(cex = 1.5),                 # X-axis title size
       ylab = list(cex = 1.5),                 # Y-axis title size
       main = list(cex = 1.8),                 # Main title size
       scales = list(x = list(cex = 1.2),      # X-axis numbers size
                     y = list(cex = 1.2)))     # Y-axis gene names size
dev.off()

library(ggplot2)

# 1. Extract and format the GBM importance data
# Assuming gbm_imp is your caret varImp object
gbm_plot_data <- data.frame(
  Gene = rownames(gbm_imp$importance),
  Importance = gbm_imp$importance$Overall
)

# Sort and take top 20
gbm_plot_data <- gbm_plot_data[order(-gbm_plot_data$Importance), ][1:20, ]
gbm_plot_data$Gene <- factor(gbm_plot_data$Gene, levels = rev(gbm_plot_data$Gene))

# 2. Create the GBM-specific plot
p_gbm <- ggplot(gbm_plot_data, aes(x = Gene, y = Importance)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#00A087FF", color = "black", size = 0.2) + # NPG Teal/Green
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12, face = "italic", color = "black"), # Larger gene names
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Top 20 Important Genes (GBM)",
    x = NULL,
    y = "Relative Variable Importance Score"
  )

# 3. Save as high-resolution TIFF
tiff("reanalyze_Figure_GBM_Importance_SciReports.tiff", width = 7, height = 8, units = "in", res = 300, compression = "lzw")
print(p_gbm)
dev.off()



library(gbm)
pls_imp <- varImp(models$PLS, scale=TRUE)
tiff("PLS_Importance.tiff", width=8, height=6, units="in", res=300)
pls_imp <- varImp(models$PLS, scale=TRUE)
# Store the plot first
p <- plot(pls_imp, top = 20, main = "Top 20 Important Genes (PLS)")

# Update with larger text sizes
update(p, 
       xlab = list(cex = 1.5),                 # X-axis title size
       ylab = list(cex = 1.5),                 # Y-axis title size
       main = list(cex = 1.8),                 # Main title size
       scales = list(x = list(cex = 1.2),      # X-axis numbers size
                     y = list(cex = 1.2)))     # Y-axis gene names size
dev.off()

library(ggplot2)

# 1. Extract and format the GBM importance data
# Assuming gbm_imp is your caret varImp object
pls_plot_data <- data.frame(
  Gene = rownames(pls_imp$importance),
  Importance = pls_imp$importance$Overall
)

# Sort and take top 20
pls_plot_data <- pls_plot_data[order(-pls_plot_data$Importance), ][1:20, ]
pls_plot_data$Gene <- factor(pls_plot_data$Gene, levels = rev(pls_plot_data$Gene))

# 2. Create the GBM-specific plot
p_pls <- ggplot(pls_plot_data, aes(x = Gene, y = Importance)) +
  geom_bar(stat = "identity", width = 0.7, fill = "lightblue", color = "black", size = 0.2) + # NPG Teal/Green
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12, face = "italic", color = "black"), # Larger gene names
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Top 20 Important Genes (PLS)",
    x = NULL,
    y = "Relative Variable Importance Score"
  )

# 3. Save as high-resolution TIFF
tiff("reanalyze_Figure_PLS_Importance.tiff", width = 7, height = 8, units = "in", res = 300, compression = "lzw")
print(p_pls)
dev.off()



# GO Enrichment
ego_up<- enrichGO(gene=selected_genes_30[1:15], universe=rownames(tt), OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont="ALL")

ego_down<- enrichGO(gene=selected_genes_30[16:30], universe=rownames(tt), OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont="ALL")

tiff("GO_up_Enrichment.tiff", width=10, height=8, units="in", res=300)
if(nrow(as.data.frame(ego_up)) > 0) print(dotplot(ego_up, showCategory=20) + ggtitle("GO: ALL (Up-regulated)"))
dev.off()

tiff("GO_down_Enrichment.tiff", width=10, height=8, units="in", res=300)
if(nrow(as.data.frame(ego_down)) > 0) print(dotplot(ego_down, showCategory=20) + ggtitle("GO: ALL (Down-regulated)"))
dev.off()

# 1. Generate Summarized Data for all Models and Datasets
library(dplyr)
library(ggplot2)


all_cm_data <- list()

for (m_name in c("PLS", "GBM")) {
  for (d_name in names(datasets)) {
    
    probs <- predict(models[[m_name]], datasets[[d_name]], type = "prob")[, "AD"]
    preds <- factor(ifelse(probs >= frozen_thresholds[[m_name]], "AD", "Control"), levels = c("AD", "Control"))
    truth <- factor(datasets[[d_name]]$condition, levels = c("AD", "Control"))
    
    # Create the counts
    cm_df <- as.data.frame(table(Prediction = preds, Target = truth))
    
    # Assign standard Confusion Matrix labels assuming "AD" is the positive class
    cm_df <- cm_df %>%
      mutate(
        Type = case_when(
          Target == "AD"      & Prediction == "AD"      ~ "TP",
          Target == "Control" & Prediction == "Control" ~ "TN",
          Target == "Control" & Prediction == "AD"      ~ "FP",
          Target == "AD"      & Prediction == "Control" ~ "FN"
        ),
        Model = m_name, 
        Dataset = d_name,
        # Dynamically grab the row count (N) for this dataset
        N_Size = nrow(datasets[[d_name]])
      )
    
    all_cm_data[[paste(m_name, d_name)]] <- cm_df
  }
}

# Combine data and append the N sizes to the header string names
plot_df <- bind_rows(all_cm_data) %>%
  mutate(
    Dataset_Label = case_when(
      Dataset == "TRAIN" ~ paste0("TRAIN (N = ", N_Size, ")"),
      Dataset == "TEST"  ~ paste0("TEST (N = ", N_Size, ")"),
      Dataset == "EXT"   ~ paste0("EXT (N = ", N_Size, ")")
    )
  )

# Extract unique labels in order to set factors correctly
unique_labels <- c(
  paste0("TRAIN (N = ", nrow(datasets$TRAIN), ")"),
  paste0("TEST (N = ", nrow(datasets$TEST), ")"),
  paste0("EXT (N = ", nrow(datasets$EXT), ")")
)

plot_df <- plot_df %>%
  mutate(Dataset_Label = factor(Dataset_Label, levels = unique_labels))

# 2. Plotting the Clean Grid with N-sizes in Headings
tiff("Figure_ConfusionMatrices_Grid.tiff", width = 10, height = 6, units = "in", res = 300)

ggplot(plot_df, aes(x = Target, y = Prediction)) +
  geom_tile(fill = "#f7f7f7", color = "grey80", linewidth = 1) +
  # Cell raw counts and types
  geom_text(aes(label = sprintf("%d\n(%s)", Freq, Type)), 
            fontface = "bold", size = 4.5, color = "black", lineheight = 0.9) +
  # Create the Grid: Facet by the new Dataset_Label containing N sizes
  facet_grid(Model ~ Dataset_Label) +
  theme_bw() +
  # Keep reference columns on top
  scale_x_discrete(position = "top") + 
  labs(title = "Confusion Matrix Grid: Model Performance Across Cohorts",
       subtitle = "Rows: Models | Columns: Datasets",
       x = "Reference (Actual)", y = "Prediction (Model)") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95"),
    axis.text = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(vjust = -0.5, face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none" 
  )

dev.off()

# =====================================================================
# 9. ENSEMBLING PLS AND GBM & UPDATING METRICS TABLE
# =====================================================================
message("Building PLS-GBM Ensemble and generating multi-dataset comparisons...")

# Define ensemble weights based on their relative Mean_CV scores
w_pls <- cv_means["PLS"] / (cv_means["PLS"] + cv_means["GBM"])
w_gbm <- cv_means["GBM"] / (cv_means["PLS"] + cv_means["GBM"])

# Extract cross-validated holdout predictions for thresholding
pred_pls <- fit_pls$pred
for(param in names(fit_pls$bestTune)) pred_pls <- pred_pls[pred_pls[[param]] == fit_pls$bestTune[[param]], ]
# Match rows by resampling fold and observation index to ensure alignment
pred_pls <- pred_pls[order(pred_pls$Resample, pred_pls$rowIndex), ]

pred_gbm <- fit_gbm$pred
for(param in names(fit_gbm$bestTune)) pred_gbm <- pred_gbm[pred_gbm[[param]] == fit_gbm$bestTune[[param]], ]
pred_gbm <- pred_gbm[order(pred_gbm$Resample, pred_gbm$rowIndex), ]

# Compute un-contaminated cross-validated ensemble probabilities
ens_cv_probs <- (w_pls * pred_pls$AD) + (w_gbm * pred_gbm$AD)
ens_roc_obj  <- pROC::roc(pred_pls$obs, ens_cv_probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
ens_threshold <- as.numeric(pROC::coords(ens_roc_obj, "best", ret="threshold", best.method="youden")[1])

# Calculate CV Mean and SD for the Ensemble
ens_cv_means_list <- sapply(unique(pred_pls$Resample), function(fold) {
  fold_idx <- which(pred_pls$Resample == fold)
  as.numeric(pROC::auc(pROC::roc(pred_pls$obs[fold_idx], ens_cv_probs[fold_idx], quiet=TRUE, levels=c("Control", "AD"), direction="<")))
})
ens_cv_mean <- mean(ens_cv_means_list)
ens_cv_sd   <- sd(ens_cv_means_list)

# Generate performance rows for the ensemble model across all datasets
ens_table_data <- list()
for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- factor(data_curr$condition, levels = c("Control", "AD"))
  
  prob_pls  <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
  prob_gbm  <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
  prob_ens  <- (w_pls * prob_pls) + (w_gbm * prob_gbm)
  
  preds_ens <- factor(ifelse(prob_ens >= ens_threshold, "AD", "Control"), levels = c("Control", "AD"))
  cm_ens    <- confusionMatrix(preds_ens, truth, positive = "AD")
  
  ens_table_data[[paste(d, "ENSEMBLE", sep = "_")]] <- data.frame(
    Model = "ENSEMBLE", Dataset = d, Accuracy = round(cm_ens$overall["Accuracy"] * 100, 1),
    Precision = round(cm_ens$byClass["Pos Pred Value"] * 100, 1), Recall = round(cm_ens$byClass["Sensitivity"] * 100, 1),
    F1_score = round(cm_ens$byClass["F1"] * 100, 1), Mean_CV = round(ens_cv_mean * 100, 1), SD_CV = round(ens_cv_sd, 3)
  )
}

# Combine with your existing table
final_output_table <- rbind(final_output_table, do.call(rbind, ens_table_data))
print(final_output_table)

# =====================================================================
# 10. GENERATING STANDALONE ROC PLOTS WITH INTEGRATED INTERNAL BOXES
# =====================================================================
message("Plotting standalone ROC curves with custom embedded legends...")

# Set up raw coordinate generation for base tracking
plot_list <- list()
for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- data_curr$condition
  
  for (m in c(names(models), "ENSEMBLE")) {
    if (m == "ENSEMBLE") {
      p_pls <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
      p_gbm <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
      probs <- (w_pls * p_pls) + (w_gbm * p_gbm)
    } else {
      probs <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    }
    
    roc_curve <- pROC::roc(truth, probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
    plot_list[[paste(d, m, sep="_")]] <- data.frame(
      Specificity = roc_curve$specificities,
      Sensitivity = roc_curve$sensitivities,
      Model = m,
      AUC = as.numeric(pROC::auc(roc_curve)),
      Dataset = d
    )
  }
}
roc_plot_df <- do.call(rbind, plot_list)

# Generate individual standalone plots with internal boxes
for (d in names(datasets)) {
  
  # Filter coordinates and output metrics tables cleanly
  cohort_df <- roc_plot_df[roc_plot_df$Dataset == d, ]
  metrics_sub <- final_output_table[final_output_table$Dataset == d, ]
  
  # Isolate and format unique AUC strings for this cohort split
  auc_lines <- c()
  unique_models <- c("GLMNET", "SVM", "PLS", "GBM", "NNET", "ENSEMBLE")
  for(m in unique_models) {
    val_auc <- round(cohort_df$AUC[cohort_df$Model == m][1], 2)
    auc_lines <- c(auc_lines, paste0(m, " (AUC = ", format(val_auc, nsmall = 2), ")"))
  }
  
  # Format the text box content
  param_box_text <- paste0(
    "--- METRICS (", d, ") ---\n",
    paste(apply(metrics_sub, 1, function(r) {
      paste0(stringr::str_pad(r["Model"], 8, "right"), " -> Acc: ", r["Accuracy"], "%, F1: ", r["F1_score"], "%")
    }), collapse = "\n"),
    "\n\n--- AREA UNDER CURVE ---\n",
    paste(auc_lines, collapse = "\n")
  )
  
  # Build plot canvas
  g_single <- ggplot(cohort_df, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
    geom_path(linewidth = 1.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = paste0("ROC Operational Performance: ", d, " Cohort"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    # Embed the custom text parameters box at the bottom-right corner
    annotate(
      "label", 
      x = 0.98, y = 0.02, 
      label = param_box_text, 
      hjust = 1, vjust = 0,
      family = "mono", size = 3.0,
      fill = "white", color = "black", alpha = 0.92,
      label.padding = unit(0.5, "lines"),
      label.size = 0.4
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "none", # Turns off standard outside legend
      panel.grid.major = element_line(color = "grey93"),
      panel.grid.minor = element_blank()
    )
  
  # Save standalone high-resolution asset file
  file_out_name <- paste0("reanalyze_Figure_ROC_Standalone_", d, ".tiff")
  tiff(file = file_out_name, unit = "in", res = 300, width = 7, height = 7, compression = "lzw")
  print(g_single)
  dev.off()
  
  message("Saved standalone asset block: ", file_out_name)
}


# =====================================================================
# 9. ENSEMBLING PLS AND GBM & UPDATING METRICS TABLE
# =====================================================================
message("Building PLS-GBM Ensemble and generating multi-dataset comparisons...")

# Define ensemble weights based on their relative Mean_CV scores
w_pls <- cv_means["PLS"] / (cv_means["PLS"] + cv_means["GBM"])
w_gbm <- cv_means["GBM"] / (cv_means["PLS"] + cv_means["GBM"])

# Extract cross-validated holdout predictions for thresholding
pred_pls <- fit_pls$pred
for(param in names(fit_pls$bestTune)) pred_pls <- pred_pls[pred_pls[[param]] == fit_pls$bestTune[[param]], ]
pred_pls <- pred_pls[order(pred_pls$Resample, pred_pls$rowIndex), ]

pred_gbm <- fit_gbm$pred
for(param in names(fit_gbm$bestTune)) pred_gbm <- pred_gbm[pred_gbm[[param]] == fit_gbm$bestTune[[param]], ]
pred_gbm <- pred_gbm[order(pred_gbm$Resample, pred_gbm$rowIndex), ]

# Compute un-contaminated cross-validated ensemble probabilities
ens_cv_probs <- (w_pls * pred_pls$AD) + (w_gbm * pred_gbm$AD)
ens_roc_obj  <- pROC::roc(pred_pls$obs, ens_cv_probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
ens_threshold <- as.numeric(pROC::coords(ens_roc_obj, "best", ret="threshold", best.method="youden")[1])

# Calculate CV Mean and SD for the Ensemble
ens_cv_means_list <- sapply(unique(pred_pls$Resample), function(fold) {
  fold_idx <- which(pred_pls$Resample == fold)
  as.numeric(pROC::auc(pROC::roc(pred_pls$obs[fold_idx], ens_cv_probs[fold_idx], quiet=TRUE, levels=c("Control", "AD"), direction="<")))
})
ens_cv_mean <- mean(ens_cv_means_list)
ens_cv_sd   <- sd(ens_cv_means_list)

# Generate performance rows for the ensemble model across all datasets
ens_table_data <- list()
for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- factor(data_curr$condition, levels = c("Control", "AD"))
  
  prob_pls  <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
  prob_gbm  <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
  prob_ens  <- (w_pls * prob_pls) + (w_gbm * prob_gbm)
  
  preds_ens <- factor(ifelse(prob_ens >= ens_threshold, "AD", "Control"), levels = c("Control", "AD"))
  cm_ens    <- confusionMatrix(preds_ens, truth, positive = "AD")
  
  ens_table_data[[paste(d, "ENSEMBLE", sep = "_")]] <- data.frame(
    Model = "ENSEMBLE", Dataset = d, Accuracy = round(cm_ens$overall["Accuracy"] * 100, 1),
    Precision = round(cm_ens$byClass["Pos Pred Value"] * 100, 1), Recall = round(cm_ens$byClass["Sensitivity"] * 100, 1),
    F1_score = round(cm_ens$byClass["F1"] * 100, 1), Mean_CV = round(ens_cv_mean * 100, 1), SD_CV = round(ens_cv_sd, 3)
  )
}

# Combine with your existing table
final_output_table <- rbind(final_output_table, do.call(rbind, ens_table_data))

# =====================================================================
# 10. GENERATING STANDALONE ROC PLOTS WITH 95% AUC CONFIDENCE INTERVALS
# =====================================================================
message("Plotting standalone ROC curves with DeLong 95% AUC confidence intervals...")

# Set up raw coordinate generation and calculate CIs using DeLong's method
plot_list <- list()
ci_list <- list()

for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- data_curr$condition
  
  for (m in c(names(models), "ENSEMBLE")) {
    if (m == "ENSEMBLE") {
      p_pls <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
      p_gbm <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
      probs <- (w_pls * p_pls) + (w_gbm * p_gbm)
    } else {
      probs <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    }
    
    roc_curve <- pROC::roc(truth, probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
    
    # Compute 95% Confidence Interval via DeLong
    roc_ci <- pROC::ci.auc(roc_curve, method="delong")
    
    plot_list[[paste(d, m, sep="_")]] <- data.frame(
      Specificity = roc_curve$specificities,
      Sensitivity = roc_curve$sensitivities,
      Model = m,
      Dataset = d
    )
    
    # Store CI limits safely
    ci_list[[paste(d, m, sep="_")]] <- data.frame(
      Dataset = d, Model = m,
      AUC = as.numeric(roc_ci[2]),
      Lower_CI = as.numeric(roc_ci[1]),
      Upper_CI = as.numeric(roc_ci[3])
    )
  }
}
roc_plot_df <- do.call(rbind, plot_list)
auc_ci_df   <- do.call(rbind, ci_list)

# Generate individual standalone plots with 95% CI parameters box
for (d in names(datasets)) {
  
  # Filter active dataset footprints
  cohort_df <- roc_plot_df[roc_plot_df$Dataset == d, ]
  metrics_sub <- final_output_table[final_output_table$Dataset == d, ]
  ci_sub      <- auc_ci_df[auc_ci_df$Dataset == d, ]
  
  # Construct formatted text lines containing Model, Accuracy, and AUC (95% CI)
  param_box_lines <- c()
  unique_models <- c("GLMNET", "SVM", "PLS", "GBM", "NNET", "ENSEMBLE")
  
  for(m in unique_models) {
    acc_val <- metrics_sub$Accuracy[metrics_sub$Model == m]
    m_ci <- ci_sub[ci_sub$Model == m, ]
    
    line_str <- paste0(
      stringr::str_pad(m, 8, "right"), 
      " Acc: ", stringr::str_pad(paste0(acc_val, "%"), 5, "right"),
      " | AUC: ", format(round(m_ci$AUC, 2), nsmall=2), 
      " (95% CI: ", format(round(m_ci$Lower_CI, 2), nsmall=2), "-", format(round(m_ci$Upper_CI, 2), nsmall=2), ")"
    )
    param_box_lines <- c(param_box_lines, line_str)
  }
  
  param_box_text <- paste0(
    "--- MODEL PERFORMANCE SUMMARY (", d, ") ---\n",
    paste(param_box_lines, collapse = "\n")
  )
  
  # Build plot canvas
  g_single <- ggplot(cohort_df, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
    geom_path(linewidth = 1.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = paste0("ROC Operational Performance: ", d, " Cohort"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    # Embed the custom confidence intervals parameters box at the bottom-right corner
    annotate(
      "label", 
      x = 0.98, y = 0.02, 
      label = param_box_text, 
      hjust = 1, vjust = 0,
      family = "mono", size = 2.8, # Slightly smaller font to fit the CI strings cleanly
      fill = "white", color = "black", alpha = 0.94,
      label.padding = unit(0.5, "lines"),
      label.size = 0.4
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "none", # Turns off standard outside legend
      panel.grid.major = element_line(color = "grey93"),
      panel.grid.minor = element_blank()
    )
  
  # Save standalone high-resolution asset file
  file_out_name <- paste0("reanalyze_Figure_ROC_Standalone_", d, ".tiff")
  tiff(file = file_out_name, unit = "in", res = 300, width = 7.5, height = 7.5, compression = "lzw")
  print(g_single)
  dev.off()
  
  message("Saved standalone asset block with Confidence Intervals: ", file_out_name)
}


# =====================================================================
# 9. ENSEMBLING PLS AND GBM & UPDATING METRICS TABLE
# =====================================================================
message("Building PLS-GBM Ensemble and generating multi-dataset comparisons...")

# Define ensemble weights based on their relative Mean_CV scores
w_pls <- cv_means["PLS"] / (cv_means["PLS"] + cv_means["GBM"])
w_gbm <- cv_means["GBM"] / (cv_means["PLS"] + cv_means["GBM"])

# Extract cross-validated holdout predictions for thresholding
pred_pls <- fit_pls$pred
for(param in names(fit_pls$bestTune)) pred_pls <- pred_pls[pred_pls[[param]] == fit_pls$bestTune[[param]], ]
pred_pls <- pred_pls[order(pred_pls$Resample, pred_pls$rowIndex), ]

pred_gbm <- fit_gbm$pred
for(param in names(fit_gbm$bestTune)) pred_gbm <- pred_gbm[pred_gbm[[param]] == fit_gbm$bestTune[[param]], ]
pred_gbm <- pred_gbm[order(pred_gbm$Resample, pred_gbm$rowIndex), ]

# Compute un-contaminated cross-validated ensemble probabilities
ens_cv_probs <- (w_pls * pred_pls$AD) + (w_gbm * pred_gbm$AD)
ens_roc_obj  <- pROC::roc(pred_pls$obs, ens_cv_probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
ens_threshold <- as.numeric(pROC::coords(ens_roc_obj, "best", ret="threshold", best.method="youden")[1])

# Calculate CV Mean and SD for the Ensemble
ens_cv_means_list <- sapply(unique(pred_pls$Resample), function(fold) {
  fold_idx <- which(pred_pls$Resample == fold)
  as.numeric(pROC::auc(pROC::roc(pred_pls$obs[fold_idx], ens_cv_probs[fold_idx], quiet=TRUE, levels=c("Control", "AD"), direction="<")))
})
ens_cv_mean <- mean(ens_cv_means_list)
ens_cv_sd   <- sd(ens_cv_means_list)

# Generate performance rows for the ensemble model across all datasets
ens_table_data <- list()
for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- factor(data_curr$condition, levels = c("Control", "AD"))
  
  prob_pls  <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
  prob_gbm  <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
  prob_ens  <- (w_pls * prob_pls) + (w_gbm * prob_gbm)
  
  preds_ens <- factor(ifelse(prob_ens >= ens_threshold, "AD", "Control"), levels = c("Control", "AD"))
  cm_ens    <- confusionMatrix(preds_ens, truth, positive = "AD")
  
  ens_table_data[[paste(d, "ENSEMBLE", sep = "_")]] <- data.frame(
    Model = "ENSEMBLE", Dataset = d, Accuracy = round(cm_ens$overall["Accuracy"] * 100, 1),
    Precision = round(cm_ens$byClass["Pos Pred Value"] * 100, 1), Recall = round(cm_ens$byClass["Sensitivity"] * 100, 1),
    F1_score = round(cm_ens$byClass["F1"] * 100, 1), Mean_CV = round(ens_cv_mean * 100, 1), SD_CV = round(ens_cv_sd, 3)
  )
}

# Combine with your existing table
final_output_table <- rbind(final_output_table, do.call(rbind, ens_table_data))

# =====================================================================
# 10. GENERATING STANDALONE ROC PLOTS WITH 95% AUC CONFIDENCE INTERVALS
# =====================================================================
message("Plotting standalone ROC curves with DeLong 95% AUC confidence intervals...")

# Set up raw coordinate generation and calculate CIs using DeLong's method
plot_list <- list()
ci_list <- list()

for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- data_curr$condition
  
  for (m in c(names(models), "ENSEMBLE")) {
    if (m == "ENSEMBLE") {
      p_pls <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
      p_gbm <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
      probs <- (w_pls * p_pls) + (w_gbm * p_gbm)
    } else {
      probs <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    }
    
    roc_curve <- pROC::roc(truth, probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
    
    # Compute 95% Confidence Interval via DeLong
    roc_ci <- pROC::ci.auc(roc_curve, method="delong")
    
    plot_list[[paste(d, m, sep="_")]] <- data.frame(
      Specificity = roc_curve$specificities,
      Sensitivity = roc_curve$sensitivities,
      Model = m,
      Dataset = d
    )
    
    # Store CI limits safely
    ci_list[[paste(d, m, sep="_")]] <- data.frame(
      Dataset = d, Model = m,
      AUC = as.numeric(roc_ci[2]),
      Lower_CI = as.numeric(roc_ci[1]),
      Upper_CI = as.numeric(roc_ci[3])
    )
  }
}
roc_plot_df <- do.call(rbind, plot_list)
auc_ci_df   <- do.call(rbind, ci_list)

# Generate individual standalone plots with 95% CI parameters box
for (d in names(datasets)) {
  
  # Filter active dataset footprints
  cohort_df <- roc_plot_df[roc_plot_df$Dataset == d, ]
  metrics_sub <- final_output_table[final_output_table$Dataset == d, ]
  ci_sub      <- auc_ci_df[auc_ci_df$Dataset == d, ]
  
  # Construct formatted text lines containing Model, Accuracy, and AUC (95% CI)
  param_box_lines <- c()
  unique_models <- c("GLMNET", "SVM", "PLS", "GBM", "NNET", "ENSEMBLE")
  
  for(m in unique_models) {
    acc_val <- metrics_sub$Accuracy[metrics_sub$Model == m]
    m_ci <- ci_sub[ci_sub$Model == m, ]
    
    line_str <- paste0(
      stringr::str_pad(m, 8, "right"), 
      " Acc: ", stringr::str_pad(paste0(acc_val, "%"), 5, "right"),
      " | AUC: ", format(round(m_ci$AUC, 2), nsmall=2), 
      " (95% CI: ", format(round(m_ci$Lower_CI, 2), nsmall=2), "-", format(round(m_ci$Upper_CI, 2), nsmall=2), ")"
    )
    param_box_lines <- c(param_box_lines, line_str)
  }
  
  param_box_text <- paste0(
    "--- MODEL PERFORMANCE SUMMARY (", d, ") ---\n",
    paste(param_box_lines, collapse = "\n")
  )
  
  # Build plot canvas
  g_single <- ggplot(cohort_df, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
    geom_path(linewidth = 1.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = paste0("ROC Operational Performance: ", d, " Cohort"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    # Embed the custom confidence intervals parameters box at the bottom-right corner
    annotate(
      "label", 
      x = 0.98, y = 0.02, 
      label = param_box_text, 
      hjust = 1, vjust = 0,
      family = "mono", size = 2.8, # Slightly smaller font to fit the CI strings cleanly
      fill = "white", color = "black", alpha = 0.94,
      label.padding = unit(0.5, "lines"),
      label.size = 0.4
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "none", # Turns off standard outside legend
      panel.grid.major = element_line(color = "grey93"),
      panel.grid.minor = element_blank()
    )
  
  # Save standalone high-resolution asset file
  file_out_name <- paste0("reanalyze_Figure_ROC_Standalone_", d, ".tiff")
  tiff(file = file_out_name, unit = "in", res = 300, width = 7.5, height = 7.5, compression = "lzw")
  print(g_single)
  dev.off()
  
  message("Saved standalone asset block with Confidence Intervals: ", file_out_name)
}


# =====================================================================
# 9. ENSEMBLING PLS AND GBM & UPDATING METRICS TABLE
# =====================================================================
message("Building PLS-GBM Ensemble and generating multi-dataset comparisons...")

# Define ensemble weights based on their relative Mean_CV scores
w_pls <- cv_means["PLS"] / (cv_means["PLS"] + cv_means["GBM"])
w_gbm <- cv_means["GBM"] / (cv_means["PLS"] + cv_means["GBM"])

# Extract cross-validated holdout predictions for thresholding
pred_pls <- fit_pls$pred
for(param in names(fit_pls$bestTune)) pred_pls <- pred_pls[pred_pls[[param]] == fit_pls$bestTune[[param]], ]
pred_pls <- pred_pls[order(pred_pls$Resample, pred_pls$rowIndex), ]

pred_gbm <- fit_gbm$pred
for(param in names(fit_gbm$bestTune)) pred_gbm <- pred_gbm[pred_gbm[[param]] == fit_gbm$bestTune[[param]], ]
pred_gbm <- pred_gbm[order(pred_gbm$Resample, pred_gbm$rowIndex), ]

# Compute un-contaminated cross-validated ensemble probabilities
ens_cv_probs <- (w_pls * pred_pls$AD) + (w_gbm * pred_gbm$AD)
ens_roc_obj  <- pROC::roc(pred_pls$obs, ens_cv_probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
ens_threshold <- as.numeric(pROC::coords(ens_roc_obj, "best", ret="threshold", best.method="youden")[1])

# Calculate CV Mean and SD for the Ensemble
ens_cv_means_list <- sapply(unique(pred_pls$Resample), function(fold) {
  fold_idx <- which(pred_pls$Resample == fold)
  as.numeric(pROC::auc(pROC::roc(pred_pls$obs[fold_idx], ens_cv_probs[fold_idx], quiet=TRUE, levels=c("Control", "AD"), direction="<")))
})
ens_cv_mean <- mean(ens_cv_means_list)
ens_cv_sd   <- sd(ens_cv_means_list)

# Generate performance rows for the ensemble model across all datasets
ens_table_data <- list()
for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- factor(data_curr$condition, levels = c("Control", "AD"))
  
  prob_pls  <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
  prob_gbm  <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
  prob_ens  <- (w_pls * prob_pls) + (w_gbm * prob_gbm)
  
  preds_ens <- factor(ifelse(prob_ens >= ens_threshold, "AD", "Control"), levels = c("Control", "AD"))
  cm_ens    <- confusionMatrix(preds_ens, truth, positive = "AD")
  
  ens_table_data[[paste(d, "ENSEMBLE", sep = "_")]] <- data.frame(
    Model = "ENSEMBLE", Dataset = d, Accuracy = round(cm_ens$overall["Accuracy"] * 100, 1),
    Precision = round(cm_ens$byClass["Pos Pred Value"] * 100, 1), Recall = round(cm_ens$byClass["Sensitivity"] * 100, 1),
    F1_score = round(cm_ens$byClass["F1"] * 100, 1), Mean_CV = round(ens_cv_mean * 100, 1), SD_CV = round(ens_cv_sd, 3)
  )
}

# Combine with your existing table
final_output_table <- rbind(final_output_table, do.call(rbind, ens_table_data))

# =====================================================================
# 10. GENERATING STANDALONE ROC PLOTS WITH 95% AUC CONFIDENCE INTERVALS
# =====================================================================
message("Plotting standalone ROC curves with DeLong 95% AUC confidence intervals...")

# Set up raw coordinate generation and calculate CIs using DeLong's method
plot_list <- list()
ci_list <- list()

for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- data_curr$condition
  
  for (m in c(names(models), "ENSEMBLE")) {
    if (m == "ENSEMBLE") {
      p_pls <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
      p_gbm <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
      probs <- (w_pls * p_pls) + (w_gbm * p_gbm)
    } else {
      probs <- predict(models[[m]], data_curr, type = "prob")[, "AD"]
    }
    
    roc_curve <- pROC::roc(truth, probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
    
    # Compute 95% Confidence Interval via DeLong
    roc_ci <- pROC::ci.auc(roc_curve, method="delong")
    
    plot_list[[paste(d, m, sep="_")]] <- data.frame(
      Specificity = roc_curve$specificities,
      Sensitivity = roc_curve$sensitivities,
      Model = m,
      Dataset = d
    )
    
    # Store CI limits safely
    ci_list[[paste(d, m, sep="_")]] <- data.frame(
      Dataset = d, Model = m,
      AUC = as.numeric(roc_ci[2]),
      Lower_CI = as.numeric(roc_ci[1]),
      Upper_CI = as.numeric(roc_ci[3])
    )
  }
}
roc_plot_df <- do.call(rbind, plot_list)
auc_ci_df   <- do.call(rbind, ci_list)

# Generate individual standalone plots with 95% CI parameters box
for (d in names(datasets)) {
  
  # Filter active dataset footprints
  cohort_df <- roc_plot_df[roc_plot_df$Dataset == d, ]
  metrics_sub <- final_output_table[final_output_table$Dataset == d, ]
  ci_sub      <- auc_ci_df[auc_ci_df$Dataset == d, ]
  
  # Construct formatted text lines containing Model, Accuracy, and AUC (95% CI)
  param_box_lines <- c()
  unique_models <- c("GLMNET", "SVM", "PLS", "GBM", "NNET", "ENSEMBLE")
  
  for(m in unique_models) {
    acc_val <- metrics_sub$Accuracy[metrics_sub$Model == m]
    m_ci <- ci_sub[ci_sub$Model == m, ]
    
    line_str <- paste0(
      stringr::str_pad(m, 8, "right"), 
      " Acc: ", stringr::str_pad(paste0(acc_val, "%"), 5, "right"),
      " | AUC: ", format(round(m_ci$AUC, 2), nsmall=2), 
      " (95% CI: ", format(round(m_ci$Lower_CI, 2), nsmall=2), "-", format(round(m_ci$Upper_CI, 2), nsmall=2), ")"
    )
    param_box_lines <- c(param_box_lines, line_str)
  }
  
  param_box_text <- paste0(
    "--- MODEL PERFORMANCE SUMMARY (", d, ") ---\n",
    paste(param_box_lines, collapse = "\n")
  )
  
  # Build plot canvas
  g_single <- ggplot(cohort_df, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
    geom_path(linewidth = 1.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = paste0("ROC Operational Performance: ", d, " Cohort"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    # Embed the custom confidence intervals parameters box at the bottom-right corner
    annotate(
      "label", 
      x = 0.98, y = 0.02, 
      label = param_box_text, 
      hjust = 1, vjust = 0,
      family = "mono", size = 2.8, # Slightly smaller font to fit the CI strings cleanly
      fill = "white", color = "black", alpha = 0.94,
      label.padding = unit(0.5, "lines"),
      label.size = 0.4
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "none", # Turns off standard outside legend
      panel.grid.major = element_line(color = "grey93"),
      panel.grid.minor = element_blank()
    )
  
  # Save standalone high-resolution asset file
  file_out_name <- paste0("reanalyze_Figure_ROC_Standalone_", d, ".tiff")
  tiff(file = file_out_name, unit = "in", res = 300, width = 7.5, height = 7.5, compression = "lzw")
  print(g_single)
  dev.off()
  
  message("Saved standalone asset block with Confidence Intervals: ", file_out_name)
}


# =====================================================================
# 11. GENERATING MULTI-PANEL FIGURE (A, B, C) FOR PLS, GBM, & ENSEMBLE
# =====================================================================
message("Generating the 3-panel publication figure (A, B, C)...")

library(grid)
library(gridExtra)

# 1. Prepare data structures to capture coordinates and DeLong CIs
panel_models <- c("PLS", "GBM", "ENSEMBLE")
plot_data_list <- list()
ci_metrics <- list()

for (d in names(datasets)) {
  data_curr <- datasets[[d]]
  truth     <- data_curr$condition
  
  for (m in panel_models) {
    if (m == "ENSEMBLE") {
      p_pls <- predict(fit_pls, data_curr, type = "prob")[, "AD"]
      p_gbm <- predict(fit_gbm, data_curr, type = "prob")[, "AD"]
      probs <- (w_pls * p_pls) + (w_gbm * p_gbm)
    } else {
      model_obj <- if(m == "PLS") fit_pls else fit_gbm
      probs <- predict(model_obj, data_curr, type = "prob")[, "AD"]
    }
    
    # Compute the ROC curve and the DeLong Confidence Intervals
    roc_curve <- pROC::roc(truth, probs, quiet=TRUE, levels=c("Control", "AD"), direction="<")
    roc_ci    <- pROC::ci.auc(roc_curve, method="delong")
    
    # Store coordinates
    plot_data_list[[paste(d, m, sep="_")]] <- data.frame(
      Specificity = roc_curve$specificities,
      Sensitivity = roc_curve$sensitivities,
      Model = m,
      Cohort = d
    )
    
    # Format the legend label string
    ci_metrics[[paste(d, m, sep="_")]] <- data.frame(
      Cohort = d, Model = m,
      Label = paste0(
        stringr::str_pad(d, 5, "right"), ": ", 
        format(round(as.numeric(roc_ci[2]), 3), nsmall = 3), " (",
        format(round(as.numeric(roc_ci[1]), 3), nsmall = 3), "-",
        format(round(as.numeric(roc_ci[3]), 3), nsmall = 3), ")"
      )
    )
  }
}

full_curve_df <- do.call(rbind, plot_data_list)
full_ci_df    <- do.call(rbind, ci_metrics)

# ----------------------------------------------------------
# Color palette
# ----------------------------------------------------------
cohort_colors <- c(
  TRAIN = "#E41A1C",
  TEST  = "#377EB8",
  EXT   = "#4DAF4A"
)



panel_plots <- list()

# ----------------------------------------------------------
# Build ROC panels
# ----------------------------------------------------------
for(i in seq_along(panel_models)){
  
  m_name <- panel_models[i]
  
  panel_letter <- c("(A)", "(B)", "(C)")[i]
  
  panel_title <- if(m_name=="ENSEMBLE")
    "Ensemble: PLS & GBM"
  else
    m_name
  
  curves_sub <- subset(full_curve_df, Model==m_name)
  labels_sub <- subset(full_ci_df, Model==m_name)
  
  # Extract specific labels to build a unified legend
  label_train <- labels_sub$Label[labels_sub$Cohort=="TRAIN"]
  label_test  <- labels_sub$Label[labels_sub$Cohort=="TEST"]
  label_ext   <- labels_sub$Label[labels_sub$Cohort=="EXT"]
  
  p <-
    ggplot(
      curves_sub,
      aes(
        x = 1-Specificity,
        y = Sensitivity,
        colour = Cohort
      )
    ) +
    
    geom_path(linewidth = 1.1) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey80"
    ) +
    
    # Combined single legend with values embedded
    scale_color_manual(
      values = cohort_colors,
      breaks = c("TRAIN", "TEST", "EXT"),
      labels = c(
        paste0("TRAIN\n(", label_train, ")"),
        paste0("TEST\n(", label_test, ")"),
        paste0("EXT\n(", label_ext, ")")
      ),
      name = "AUC (95% CI)"
    ) +
    
    scale_x_continuous(
      limits = c(-0.01,1.01),
      breaks = seq(0,1,0.2)
    ) +
    
    scale_y_continuous(
      limits = c(-0.01,1.01),
      breaks = seq(0,1,0.2)
    ) +
    
    labs(
      title = panel_title,
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    
    theme_bw() +
    
    theme(
      
      plot.title =
        element_text(
          face="bold",
          size=15,
          hjust=.5,
          margin=margin(b=8)
        ),
      
      panel.grid =
        element_blank(),
      
      panel.border =
        element_rect(
          colour="black",
          fill=NA,
          linewidth=.8
        ),
      
      axis.title =
        element_text(size=12),
      
      axis.text =
        element_text(
          size=11,
          colour="black"
        ),
      
      # Adjusted single legend placement inside plot
      legend.position = c(0.52, 0.04),
      legend.justification = c(0, 0),
      
      legend.background =
        element_rect(
          fill=alpha("white", 0.95),
          colour="black",
          linewidth=0.4
        ),
      
      legend.key = element_blank(),
      
      # Smaller, compact font settings for the legend elements
      legend.title = element_text(size = 8.5, face = "bold"),
      legend.text  = element_text(size = 7.5, face = "bold"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.margin = margin(t = 4, r = 6, b = 4, l = 6)
    ) +
    # Force lines in the legend to display keys cleanly with tighter spacing
    guides(color = guide_legend(byrow = TRUE, keyheight = unit(0.5, "cm")))
  
  panel_plots[[m_name]] <-
    arrangeGrob(
      p,
      top=textGrob(
        panel_letter,
        x=unit(.05,"npc"),
        y=unit(.82,"npc"),
        just=c("left","top"),
        gp=gpar(
          fontsize=18,
          fontface="bold"
        )
      )
    )
}

# ----------------------------------------------------------
# Save figure
# ----------------------------------------------------------
tiff(
  "reanalyze_Figure_Selected_Models_Panels.tiff",
  width = 13.5,
  height = 4.8,
  units = "in",
  res = 300,
  compression = "lzw"
)

grid.arrange(
  grobs = panel_plots,
  ncol = 3
)

dev.off()
message("Saved custom isolated panel figure successfully.")

########CONFUSİON MATRİCES PLOT

##------------------------------------------------------------
## Function to compute confusion matrix counts
##------------------------------------------------------------

get_cm_data <- function(model, data, threshold, model_name, dataset_name){
  
  ## predicted probability
  prob <- predict(model, data, type = "prob")[,"AD"]
  
  ## predicted class
  pred <- factor(
    ifelse(prob >= threshold, "AD", "Control"),
    levels = c("AD","Control")
  )
  
  ## true class
  truth <- factor(
    data$condition,
    levels = c("AD","Control")
  )
  
  ## confusion matrix
  cm <- confusionMatrix(
    pred,
    truth,
    positive = "AD"
  )
  
  ## extract counts
  tbl <- as.data.frame(cm$table)
  
  names(tbl) <- c("Prediction","Target","Freq")
  
  tbl <- tbl %>%
    mutate(
      Type = case_when(
        Prediction == "AD"      & Target == "AD"      ~ "TP",
        Prediction == "Control" & Target == "Control" ~ "TN",
        Prediction == "AD"      & Target == "Control" ~ "FP",
        Prediction == "Control" & Target == "AD"      ~ "FN"
      ),
      Model   = model_name,
      Dataset = dataset_name,
      N = nrow(data)
    )
  
  return(tbl)
}

models_to_plot <- c("PLS","GBM")

plot_df <- map_dfr(models_to_plot, function(m){
  
  map_dfr(names(datasets), function(d){
    
    get_cm_data(
      model        = models[[m]],
      data         = datasets[[d]],
      threshold    = frozen_thresholds[[m]],
      model_name   = m,
      dataset_name = d
    )
    
  })
  
})

plot_df <- plot_df %>%
  mutate(
    Dataset_Label = factor(
      paste0(Dataset," (N = ",N,")"),
      levels = c(
        paste0("TRAIN (N = ",nrow(datasets$TRAIN),")"),
        paste0("TEST (N = ",nrow(datasets$TEST),")"),
        paste0("EXT (N = ",nrow(datasets$EXT),")")
      )
    )
  )

plot_df %>%
  select(Model, Dataset, Target, Prediction, Freq) %>%
  arrange(Model, Dataset)


plot_df %>%
  group_by(Model, Dataset, Target) %>%
  summarise(
    Total = sum(Freq),
    .groups = "drop"
  )

tiff("Figure_ConfusionMatrices_Grid.tiff", width = 8, height = 10, units = "in", res = 300)
ggplot(plot_df,
       aes(Target, Prediction)) +
  
  geom_tile(fill = "grey97",
            colour = "grey70",
            linewidth = 0.8) +
  
  geom_text(aes(label = paste0(Freq, "\n(", Type, ")")),
            fontface = "bold",
            size = 5.5,
            lineheight = 0.9) +
  
  facet_grid(Model ~ Dataset_Label) +
  
  scale_x_discrete(position = "top") +
  
  labs(
    x = "Reference (Actual)",
    y = "Prediction (Model)"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    # Facet labels (PLS, GBM, TRAIN, TEST, EXT)
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    
    # Axis titles
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(b = 12)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 12)),
    
    # AD / Control labels
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    
    # Remove legend
    legend.position = "none",
    
    # Panel spacing
    panel.spacing = unit(1.0, "lines")
  )
dev.off()


#=========================================================
# Function to calculate performance metrics
#=========================================================

calc_metrics <- function(model, data, threshold, model_name, dataset_name){
  
  ## Predicted probabilities
  probs <- predict(model, data, type = "prob")[,"AD"]
  
  ## Predicted classes
  pred <- factor(
    ifelse(probs >= threshold, "AD", "Control"),
    levels = c("AD","Control")
  )
  
  ## True classes
  truth <- factor(
    data$condition,
    levels = c("AD","Control")
  )
  
  ## Confusion matrix
  cm <- confusionMatrix(
    pred,
    truth,
    positive = "AD"
  )
  
  tibble(
    Dataset = dataset_name,
    Model   = model_name,
    Accuracy  = cm$overall["Accuracy"] * 100,
    Precision = cm$byClass["Pos Pred Value"] * 100,
    Recall    = cm$byClass["Sensitivity"] * 100,
    F1_score  = (2 *
                   cm$byClass["Pos Pred Value"] *
                   cm$byClass["Sensitivity"]) /
      (cm$byClass["Pos Pred Value"] +
         cm$byClass["Sensitivity"]) * 100
  )
}

#=========================================================
# Calculate metrics for all models and datasets
#=========================================================

metrics_df <-
  map_dfr(names(models), function(m){
    
    map_dfr(names(datasets), function(d){
      
      calc_metrics(
        model        = models[[m]],
        data         = datasets[[d]],
        threshold    = frozen_thresholds[[m]],
        model_name   = m,
        dataset_name = d
      )
      
    })
    
  })

#=========================================================
# Convert to long format
#=========================================================

plot_df <-
  metrics_df %>%
  pivot_longer(
    cols = c(Accuracy, Precision, Recall, F1_score),
    names_to = "Metric",
    values_to = "Score"
  )

# Order factors
plot_df$Dataset <-
  factor(plot_df$Dataset,
         levels = c("TRAIN","TEST","EXT"))

plot_df$Metric <-
  factor(plot_df$Metric,
         levels = c("Accuracy",
                    "Precision",
                    "Recall",
                    "F1_score"))

plot_df$Model <-
  factor(plot_df$Model,
         levels = c("GBM",
                    "GLMNET",
                    "NNET",
                    "PLS",
                    "SVM"))

# Colors
model_colors <- c(
  GBM     = "#4DBBD5",
  GLMNET  = "#E64B35",
  NNET    = "#7E62B3",
  PLS     = "#FDB462",
  SVM     = "#2E9F55"
)

#=========================================================
# Plot
#=========================================================

ggplot(plot_df,
       aes(Metric,
           Score,
           fill = Model)) +
  
  geom_col(position = position_dodge(0.8),
           width = 0.75) +
  
  facet_grid(Dataset ~ .) +
  
  scale_fill_manual(values = model_colors) +
  
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100,20),
    expand = expansion(mult = c(0,0.02))
  ) +
  
  labs(
    x = "Model Metrics",
    y = "Score (%)",
    fill = "Models"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    strip.text = element_text(size = 15,
                              face = "bold"),
    
    strip.background = element_rect(fill = "grey95"),
    
    axis.title = element_text(size = 16,
                              face = "bold"),
    
    axis.text.x = element_text(size = 13,
                               face = "bold"),
    
    axis.text.y = element_text(size = 12,
                               face = "bold"),
    
    legend.title = element_text(size = 13,
                                face = "bold"),
    
    legend.text = element_text(size = 12),
    
    panel.grid.major.x = element_blank(),
    
    legend.position = "bottom"
  )

# Save
ggsave(
  "Figure_ModelPerformance_Vertical_Grid.tiff",
  width = 8,
  height = 11,
  dpi = 600,
  compression = "lzw"
)

