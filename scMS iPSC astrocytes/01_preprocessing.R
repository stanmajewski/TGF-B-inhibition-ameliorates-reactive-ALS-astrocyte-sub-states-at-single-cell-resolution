#!/usr/bin/env Rscript

# =============================================================================
# SINGLE CELL PROTEOMICS PREPROCESSING PIPELINE
# =============================================================================
# This script performs quality control, filtering, and normalization of
# single-cell mass spectrometry data
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "data/raw/proteomics_data.xlsx",
              help = "Path to raw proteomics Excel file [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "results/preprocessing",
              help = "Output directory [default: %default]"),
  make_option("--min-protein-detection-global", type = "double", default = 0.10,
              help = "Minimum global protein detection rate [default: %default]"),
  make_option("--min-protein-detection-per-group", type = "double", default = 0.20,
              help = "Minimum protein detection rate per group [default: %default]"),
  make_option("--min-protein-detection-per-batch", type = "double", default = 0.10,
              help = "Minimum protein detection rate per batch [default: %default]"),
  make_option("--min-batches-detected", type = "integer", default = 3,
              help = "Minimum number of batches where protein must be detected [default: %default]"),
  make_option("--max-protein-missing", type = "double", default = 0.65,
              help = "Maximum allowed missing rate for proteins [default: %default]"),
  make_option("--cv-threshold", type = "double", default = 3.0,
              help = "Maximum coefficient of variation threshold [default: %default]"),
  make_option("--min-cell-detection", type = "double", default = 0.20,
              help = "Minimum detection rate for cell filtering [default: %default]"),
  make_option("--min-cells-per-group", type = "integer", default = 5,
              help = "Minimum cells required per group [default: %default]"),
  make_option("--log-base", type = "integer", default = 2,
              help = "Base for log transformation [default: %default]"),
  make_option("--pseudocount", type = "double", default = 0.01,
              help = "Pseudocount for log transformation [default: %default]"),
  make_option("--imputation-method", type = "character", default = "hybrid",
              help = "Imputation method (hybrid or knn) [default: %default]"),
  make_option("--knn-k", type = "integer", default = 5,
              help = "Number of neighbors for KNN imputation [default: %default]"),
  make_option("--save-checkpoints", action = "store_true", default = TRUE,
              help = "Save intermediate checkpoints [default: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nSingle-cell proteomics preprocessing pipeline",
  epilogue = "Example: Rscript 01_preprocessing.R -i data/raw/data.xlsx -o results/preprocessing"
)

opt <- parse_args(opt_parser)

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
  data_file = opt$input,
  output_dir = opt$output,
  
  filtering = list(
    min_protein_detection_global = opt$`min-protein-detection-global`,
    min_protein_detection_per_group = opt$`min-protein-detection-per-group`,
    min_protein_detection_per_batch = opt$`min-protein-detection-per-batch`,
    min_batches_detected = opt$`min-batches-detected`,
    max_protein_missing_overall = opt$`max-protein-missing`,
    cv_threshold = opt$`cv-threshold`,
    min_cell_detection = opt$`min-cell-detection`,
    min_cells_per_group = opt$`min-cells-per-group`
  ),
  
  transformation = list(
    log_base = opt$`log-base`,
    pseudocount = opt$pseudocount
  ),
  
  imputation = list(
    method = opt$`imputation-method`,
    impute_by = "Batch",
    low_missing_threshold = 0.20,
    medium_missing_threshold = 0.50,
    knn_k = opt$`knn-k`
  ),
  
  normalization = list(
    method = "median_center"
  ),
  
  save_checkpoints = opt$`save-checkpoints`
)

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("SINGLE CELL PROTEOMICS PREPROCESSING\n")
cat("=============================================================================\n\n")

cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(impute)
})

# Create output directories
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "checkpoints"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

cat("✓ Packages loaded\n\n")

# Print configuration
cat("Configuration:\n")
cat("  Input file:", CONFIG$data_file, "\n")
cat("  Output directory:", CONFIG$output_dir, "\n")
cat("  Save checkpoints:", CONFIG$save_checkpoints, "\n\n")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

save_checkpoint <- function(object, name) {
  if(CONFIG$save_checkpoints) {
    checkpoint_file <- file.path(CONFIG$output_dir, "checkpoints", 
                                 paste0(name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
    saveRDS(object, checkpoint_file)
    cat("  💾 Checkpoint saved:", basename(checkpoint_file), "\n")
  }
}

calculate_cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

# =============================================================================
# STEP 1: DATA LOADING
# =============================================================================

cat("=== STEP 1: DATA LOADING ===\n")

raw_data_temp <- read_excel(CONFIG$data_file, col_names = TRUE) %>% 
  as.data.frame()

protein_ids <- raw_data_temp[, 1]
intensity_data <- raw_data_temp[, 6:ncol(raw_data_temp)]

# Clean cell names
cell_names <- colnames(intensity_data)
cell_names <- gsub("^\\[\\d+\\]\\s*", "", cell_names)
cell_names <- gsub("\\.d(\\..*)?$", "", cell_names)  # e.g. "CTRL1_U_cell01.d.PG.Quantity" -> "CTRL1_U_cell01"
colnames(intensity_data) <- cell_names

# Clean protein IDs
protein_ids <- as.character(protein_ids)
problematic <- is.na(protein_ids) | protein_ids == "" | protein_ids == "NA"
if(any(problematic)) {
  protein_ids[problematic] <- paste0("Unknown_Protein_", which(problematic))
}
protein_ids <- make.unique(protein_ids, sep = "_")

# Create data matrix
raw_data <- as.matrix(intensity_data)
raw_data[raw_data == "NaN"] <- NA
raw_data <- apply(raw_data, 2, as.numeric)
rownames(raw_data) <- protein_ids

cat("✓ Raw data loaded:", nrow(raw_data), "proteins ×", ncol(raw_data), "cells\n\n")

# =============================================================================
# STEP 2: METADATA CREATION
# =============================================================================

cat("=== STEP 2: METADATA CREATION ===\n")

# Metadata is parsed directly from the cell names, which follow the format
#   <Line>_<Treatment>_cell<NN>   e.g. CTRL1_U_cell01, Mut1-1_SB_cell12
# Treatment: U = untreated (vehicle), SB = SB431542
# Genotype: lines starting "Mut" are VCP mutant, all others are control

name_pattern <- "^(.+)_(U|SB)_cell(\\d+)$"
bad_names <- cell_names[!grepl(name_pattern, cell_names)]
if(length(bad_names) > 0) {
  stop("Cell names not in <Line>_<U|SB>_cell<NN> format: ",
       paste(head(bad_names), collapse = ", "))
}

metadata <- data.frame(
  Cell_ID = cell_names,
  Line = sub(name_pattern, "\\1", cell_names),
  Treatment = sub(name_pattern, "\\2", cell_names),
  Cell_Number = as.integer(sub(name_pattern, "\\3", cell_names)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Genotype = ifelse(grepl("^Mut", Line), "VCP", "CTRL"),
    Group = paste0(Genotype, "_", Treatment),
    Line_Treatment = paste0(Line, "_", Treatment),
    Batch = Line   # plate info is no longer in the names; each line is used as the batch
  )

cat("✓ Metadata created for", nrow(metadata), "cells\n")
print(table(metadata$Line, metadata$Treatment))
cat("\n")

# =============================================================================
# STEP 3: LOG TRANSFORMATION
# =============================================================================

cat("=== STEP 3: LOG TRANSFORMATION ===\n")

raw_data[raw_data == 0] <- NA
log_data <- log2(raw_data + CONFIG$transformation$pseudocount)

cat("✓ Log2 transformation complete\n\n")

# =============================================================================
# STEP 4: CONTAMINANT REMOVAL
# =============================================================================

cat("=== STEP 4: CONTAMINANT REMOVAL ===\n")

protein_names <- rownames(log_data)

contaminant_patterns <- c(
  "SWISS-PROT", "SWISS\\.PROT", "TREMBL",
  ":P[0-9]{5}", ":Q[0-9]{5}",
  "^CONTAM", "^REV_", "^Reverse",
  ";KRT[0-9]", "^KRT[0-9].*;",
  "Trypsin"
)

is_contaminant <- rep(FALSE, length(protein_names))

for(pattern in contaminant_patterns) {
  matches <- grepl(pattern, protein_names, ignore.case = TRUE)
  if(any(matches)) {
    cat("  Found", sum(matches), "entries matching:", pattern, "\n")
    is_contaminant <- is_contaminant | matches
  }
}

n_contaminants <- sum(is_contaminant)
log_data <- log_data[!is_contaminant, ]

cat("✓ Removed", n_contaminants, "contaminant entries\n")
cat("  Remaining proteins:", nrow(log_data), "\n\n")

save_checkpoint(log_data, "01_contaminants_removed")

# =============================================================================
# STEP 5: QUALITY FILTERING
# =============================================================================

cat("=== STEP 5: QUALITY FILTERING ===\n")

# Protein detection rates
protein_detection_global <- rowSums(!is.na(log_data)) / ncol(log_data)
cell_detection <- colSums(!is.na(log_data)) / nrow(log_data)

# Detection by group
protein_detection_by_group <- matrix(NA, nrow = nrow(log_data), 
                                     ncol = length(unique(metadata$Group)))
colnames(protein_detection_by_group) <- unique(metadata$Group)
rownames(protein_detection_by_group) <- rownames(log_data)

for(grp in unique(metadata$Group)) {
  grp_cells <- metadata$Group == grp
  protein_detection_by_group[, grp] <- rowSums(!is.na(log_data[, grp_cells])) / sum(grp_cells)
}

# Detection by batch
protein_detection_by_batch <- matrix(NA, nrow = nrow(log_data),
                                     ncol = length(unique(metadata$Batch)))
colnames(protein_detection_by_batch) <- unique(metadata$Batch)
rownames(protein_detection_by_batch) <- rownames(log_data)

for(batch in unique(metadata$Batch)) {
  batch_cells <- metadata$Batch == batch
  protein_detection_by_batch[, batch] <- rowSums(!is.na(log_data[, batch_cells])) / sum(batch_cells)
}

# Apply filtering criteria
pass_global <- protein_detection_global >= CONFIG$filtering$min_protein_detection_global
pass_group <- apply(protein_detection_by_group, 1, max) >= CONFIG$filtering$min_protein_detection_per_group
batch_pass_threshold <- protein_detection_by_batch >= CONFIG$filtering$min_protein_detection_per_batch
batches_passed <- rowSums(batch_pass_threshold)
pass_batch <- batches_passed >= CONFIG$filtering$min_batches_detected
protein_missing_rate <- rowSums(is.na(log_data)) / ncol(log_data)
pass_missing <- protein_missing_rate <= CONFIG$filtering$max_protein_missing_overall
protein_cv <- apply(log_data, 1, calculate_cv)
pass_cv <- protein_cv < CONFIG$filtering$cv_threshold

proteins_final <- pass_global & pass_group & pass_batch & pass_missing & pass_cv
proteins_final[is.na(proteins_final)] <- FALSE

# Cell filtering
cells_to_keep <- cell_detection >= CONFIG$filtering$min_cell_detection
group_sizes_after <- table(metadata$Group[cells_to_keep])
small_groups <- names(group_sizes_after)[group_sizes_after < CONFIG$filtering$min_cells_per_group]
if(length(small_groups) > 0) {
  cells_to_keep <- cells_to_keep & !(metadata$Group %in% small_groups)
}

filtered_data <- log_data[proteins_final, cells_to_keep]
filtered_metadata <- metadata[cells_to_keep, ]

cat("✓ Proteins retained:", sum(proteins_final), "/", length(proteins_final), "\n")
cat("✓ Cells retained:", sum(cells_to_keep), "/", length(cells_to_keep), "\n\n")

save_checkpoint(list(data = filtered_data, metadata = filtered_metadata), "02_filtered")

# =============================================================================
# STEP 6: IMPUTATION
# =============================================================================

cat("=== STEP 6: HYBRID IMPUTATION ===\n")

impute_hybrid_batch_aware <- function(data_matrix, metadata_df, batch_column) {
  batches <- unique(metadata_df[[batch_column]])
  imputed_data <- data_matrix
  
  for(batch in batches) {
    cat("  Processing", batch, "...\n")
    batch_mask <- metadata_df[[batch_column]] == batch
    batch_data <- data_matrix[, batch_mask, drop = FALSE]
    protein_missing_rate <- rowSums(is.na(batch_data)) / ncol(batch_data)
    
    # Low missing: KNN imputation
    low_missing <- protein_missing_rate < CONFIG$imputation$low_missing_threshold
    med_missing <- protein_missing_rate >= CONFIG$imputation$low_missing_threshold & 
      protein_missing_rate <= CONFIG$imputation$medium_missing_threshold
    
    if(sum(low_missing) > 0) {
      low_data <- batch_data[low_missing, , drop = FALSE]
      if(any(is.na(low_data))) {
        tryCatch({
          imputed_low <- impute::impute.knn(low_data, k = CONFIG$imputation$knn_k)$data
          imputed_data[low_missing, batch_mask] <- imputed_low
        }, error = function(e) {
          for(i in which(low_missing)) {
            na_mask <- is.na(imputed_data[i, batch_mask])
            if(any(na_mask)) {
              imputed_data[i, batch_mask][na_mask] <- mean(batch_data[i, ], na.rm = TRUE)
            }
          }
        })
      }
    }
    
    # Medium missing: Left-censored imputation
    if(sum(med_missing) > 0) {
      for(i in which(med_missing)) {
        na_mask <- is.na(imputed_data[i, batch_mask])
        if(any(na_mask) && sum(!na_mask) >= 2) {
          observed <- batch_data[i, !na_mask]
          min_observed <- min(observed, na.rm = TRUE)
          obs_sd <- sd(observed, na.rm = TRUE)
          noise <- rnorm(sum(na_mask), mean = -0.5 * obs_sd, sd = 0.3 * obs_sd)
          imputed_values <- pmax(min_observed + noise, min_observed - 2 * obs_sd)
          imputed_data[i, batch_mask][na_mask] <- imputed_values
        }
      }
    }
  }
  
  # Final cleanup: impute any remaining NAs with 5th percentile
  for(i in 1:nrow(imputed_data)) {
    na_mask <- is.na(imputed_data[i, ])
    if(any(na_mask)) {
      observed <- imputed_data[i, !na_mask]
      if(length(observed) > 0) {
        fill_value <- quantile(observed, 0.05, na.rm = TRUE)
        imputed_data[i, na_mask] <- fill_value
      }
    }
  }
  
  return(imputed_data)
}

imputed_data <- impute_hybrid_batch_aware(filtered_data, filtered_metadata, CONFIG$imputation$impute_by)

cat("✓ Imputation complete\n")
cat("  Remaining NAs:", sum(is.na(imputed_data)), "\n\n")

save_checkpoint(imputed_data, "03_imputed")

# =============================================================================
# STEP 7: MEDIAN-CENTER NORMALIZATION
# =============================================================================

cat("=== STEP 7: NORMALIZATION ===\n")

col_medians <- apply(imputed_data, 2, median)
global_median <- median(col_medians)
normalized_data <- sweep(imputed_data, 2, col_medians - global_median, FUN = "-")

cat("✓ Median-center normalization complete\n\n")

save_checkpoint(normalized_data, "04_normalized")

# =============================================================================
# STEP 8: FINAL EXPORT
# =============================================================================

cat("=== STEP 8: FINAL EXPORT ===\n")

seurat_data <- normalized_data
seurat_metadata <- as.data.frame(filtered_metadata)
rownames(seurat_metadata) <- seurat_metadata$Cell_ID

final_output <- list(
  seurat_data = seurat_data,
  seurat_metadata = seurat_metadata,
  config = CONFIG,
  n_contaminants_removed = n_contaminants
)

save_checkpoint(final_output, "05_final_export")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("PREPROCESSING COMPLETE\n")
cat("=============================================================================\n\n")
cat("Summary:\n")
cat("  Input proteins:", nrow(raw_data), "\n")
cat("  Contaminants removed:", n_contaminants, "\n")
cat("  Final proteins:", nrow(seurat_data), "\n")
cat("  Final cells:", ncol(seurat_data), "\n")
cat("  Missing values:", sum(is.na(seurat_data)), "\n\n")
cat("Output directory:", CONFIG$output_dir, "\n")
cat("  └── checkpoints/05_final_export_*.rds\n")
cat("  └── tables/\n\n")
cat("Next step: Run 02_harmony_clustering.R\n\n")
