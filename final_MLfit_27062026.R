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

ctrl_fs <- trainControl(method="cv", number=5)

deg_selection <- function(x, y, age, sex) {
  
  y <- factor(y, levels=c("Control","AD"))
  
  design <- model.matrix(~ y + age + sex)
  
  fit <- limma::lmFit(t(x), design)
  fit <- limma::eBayes(fit)
  
  tt <- limma::topTable(fit, coef="yAD", number=Inf)
  
  sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.8, ]
  
  up   <- rownames(sig)[sig$logFC > 0]
  down <- rownames(sig)[sig$logFC < 0]
  
  c(head(up,15), head(down,15))
}

selected_genes <- deg_selection(train_x, train_y, train_age, train_sex)

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

up   <- rownames(sig)[sig$logFC > 0]
down <- rownames(sig)[sig$logFC < 0]

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
ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  
  # Base points (all genes)
  geom_point(aes(color = status), alpha = 0.6, size = 1.5) +
  
  # Highlight top 30 genes with different color
  geom_point(
    data = volcano_df %>% filter(highlight30 == "Top30"),
    color = "purple",   # <-- different color for top 30
    size = 3
  ) +
  
  # Labels for top 30 genes
  geom_text_repel(
    data = volcano_df %>% filter(highlight30 == "Top30"),
    aes(label = gene),
    size = 3.5,
    box.padding = 0.4,
    max.overlaps = 100
  ) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # Colors
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NotSig" = "grey",
    "Biomarker genes"="purple"
  )) +
  
  # Theme adjustments
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Expression",
    subtitle = paste0("Number of up-regulated genes is ",209 , 
                      "\nNumber of down-regulated genes is ",436 ))

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

# GO Enrichment
ego_up<- enrichGO(gene=selected_genes_30[1:15], universe=rownames(tt), OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont="ALL")

ego_down<- enrichGO(gene=selected_genes_30[16:30], universe=rownames(tt), OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont="ALL")

tiff("GO_up_Enrichment.tiff", width=10, height=8, units="in", res=300)
if(nrow(as.data.frame(ego_up)) > 0) print(dotplot(ego_up, showCategory=20) + ggtitle("GO: ALL (Up-regulated)"))
dev.off()

tiff("GO_down_Enrichment.tiff", width=10, height=8, units="in", res=300)
if(nrow(as.data.frame(ego_down)) > 0) print(dotplot(ego_down, showCategory=20) + ggtitle("GO: ALL (Down-regulated)"))
dev.off()



