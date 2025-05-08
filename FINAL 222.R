# ================================
# STT 303 Final Project: GTEx Gene Expression Analysis
# Optimized & Memory-Safe Pipeline
# ================================

# --- Load Libraries ---
library(tidyverse)
library(readr)
library(caret)
library(pROC)
library(cluster)
library(factoextra)
library(heatmaply)
library(gridExtra)
library(matrixStats)

# --- 1. Function to Read First N Genes Efficiently ---
read_gene_subset <- function(path, n_genes = 10000) {
  header_line <- read_lines(path, skip = 1, n_max = 1)
  col_names <- str_split(header_line, "\t")[[1]]
  df <- read_tsv(path, skip = 2, n_max = n_genes, col_names = col_names, show_col_types = FALSE)
  return(df)
}

# --- 2. Load Data ---
gene_expr <- read_gene_subset("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", n_genes = 10000)
metadata <- read_tsv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", show_col_types = FALSE)

# --- 3. Prepare Expression Matrix ---
first_col <- colnames(gene_expr)[1]  # Should be "Name"
gene_matrix <- gene_expr %>%
  column_to_rownames(var = first_col)

# --- 4. Ensure Numeric and Filter Low-Expression Genes ---
gene_matrix[] <- lapply(gene_matrix, function(x) as.numeric(as.character(x)))
gene_matrix <- gene_matrix[rowMedians(as.matrix(gene_matrix), na.rm = TRUE) > 1, ]

# --- 5. Transpose and Merge ---
expr_t <- as.data.frame(t(gene_matrix))
expr_t$SAMPID <- rownames(expr_t)
merged_data <- merge(expr_t, metadata, by = "SAMPID")

# --- 6. Log-transform Gene Expression ---
gene_cols <- merged_data %>%
  select(where(is.numeric)) %>%
  colnames()
merged_data[gene_cols] <- log1p(merged_data[gene_cols])

# --- 7. PCA (if enough variation) ---
pca_data <- merged_data[gene_cols]
non_constant_cols <- apply(pca_data, 2, function(x) var(x, na.rm = TRUE)) != 0

if (sum(non_constant_cols, na.rm = TRUE) >= 2) {
  pca_data <- pca_data[, non_constant_cols]
  pca_result <- prcomp(pca_data, scale. = TRUE)
  fviz_pca_ind(pca_result,
               col.ind = merged_data$SMTS,
               addEllipses = TRUE,
               label = "none") +
    ggtitle("PCA of Gene Expression (Filtered)")
} else {
  message("Not enough variable genes for PCA — try increasing n_genes.")
}

# --- 8. Heatmap of First 30 Genes (Optional) ---
if (length(gene_cols) >= 30) {
  heatmap_data <- merged_data[gene_cols[1:30]]
  heatmap_data <- heatmap_data[, colSums(is.na(heatmap_data)) == 0]
  
  if (ncol(heatmap_data) >= 2 && nrow(heatmap_data) >= 2) {
    heatmaply(heatmap_data,
              k_row = 3,
              k_col = 2,
              main = "Heatmap of First 30 Genes")
  } else {
    message("Not enough valid data to plot heatmap.")
  }
} else {
  message("Fewer than 30 gene columns — skipping heatmap.")
}

# --- 9. Optional: Regression on AGE (if available) ---
if ("AGE" %in% colnames(merged_data)) {
  merged_data$AGE <- as.numeric(factor(merged_data$AGE))
  if (length(unique(merged_data$AGE)) > 1) {
    set.seed(123)
    train_idx <- createDataPartition(merged_data$AGE, p = 0.8, list = FALSE)
    train_data <- merged_data[train_idx, ]
    test_data <- merged_data[-train_idx, ]
    model <- train(AGE ~ ., data = train_data[, c("AGE", gene_cols)], method = "rf")
    predictions <- predict(model, test_data)
    print(postResample(predictions, test_data$AGE))
  } else {
    message("AGE has only one class — skipping regression.")
  }
}

