#!/usr/bin/env Rscript

# =============================================================================
# HARMONY INTEGRATION & CLUSTERING
# =============================================================================
# This script performs batch correction using Harmony, dimensionality reduction,
# and clustering of single-cell proteomics data
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type = "character", 
              default = "results/preprocessing/checkpoints",
              help = "Path to preprocessing checkpoint directory [default: %default]"),
  make_option(c("-o", "--output"), type = "character", 
              default = "results/harmony",
              help = "Output directory [default: %default]"),
  make_option("--batch-variable", type = "character", default = "Line_Treatment",
              help = "Metadata column for batch correction [default: %default]"),
  make_option("--theta", type = "double", default = 1,
              help = "Harmony diversity clustering penalty [default: %default]"),
  make_option("--sigma", type = "double", default = 0.3,
              help = "Harmony ridge regression penalty [default: %default]"),
  make_option("--lambda", type = "double", default = 0.5,
              help = "Harmony entropy regularization [default: %default]"),
  make_option("--n-pcs", type = "integer", default = 30,
              help = "Number of PCs to use [default: %default]"),
  make_option("--resolution", type = "double", default = 0.3,
              help = "Clustering resolution [default: %default]"),
  make_option("--n-neighbors", type = "integer", default = 50,
              help = "Number of neighbors for UMAP [default: %default]"),
  make_option("--min-dist", type = "double", default = 0.5,
              help = "UMAP minimum distance [default: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nHarmony batch correction and clustering pipeline",
  epilogue = "Example: Rscript 02_harmony_clustering.R -i results/preprocessing/checkpoints -o results/harmony"
)

opt <- parse_args(opt_parser)

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
  preprocessing_checkpoint = opt$input,
  output_dir = opt$output,
  
  batch_variable = opt$`batch-variable`,
  theta = opt$theta,
  sigma = opt$sigma,
  lambda = opt$lambda,
  nclust = 50,
  max_iter = 10,
  n_pcs = opt$`n-pcs`,
  
  target_resolution = opt$resolution,
  
  n_neighbors = opt$`n-neighbors`,
  min_dist = opt$`min-dist`,
  spread = 1.5,
  metric = "euclidean"
)

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("HARMONY INTEGRATION & CLUSTERING\n")
cat("=============================================================================\n\n")

cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(colorspace)
})

# Create output directories
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

cat("✓ Packages loaded\n\n")

# Print configuration
cat("Configuration:\n")
cat("  Input checkpoint:", CONFIG$preprocessing_checkpoint, "\n")
cat("  Output directory:", CONFIG$output_dir, "\n")
cat("  Batch variable:", CONFIG$batch_variable, "\n")
cat("  Harmony theta:", CONFIG$theta, "\n")
cat("  Resolution:", CONFIG$target_resolution, "\n\n")

# Custom palette function
pal <- function(v) {
  v2 <- length(unique(v))
  if (v2 == 2) {
    p2 <- c("grey", "blue")
  } else if (v2 == 3) {
    p2 <- c("blue", "grey", "orange")
  } else if (v2 < 6) {
    p2 <- hcl.colors(6, "Set 2")[1:v2]
  } else {
    p2 <- hcl.colors(v2, "Set 2")
  }
  p2 <- desaturate(p2, amount = 0.3)
  return(p2)
}

# =============================================================================
# STEP 1: LOAD PREPROCESSED DATA
# =============================================================================

cat("=== STEP 1: LOADING PREPROCESSED DATA ===\n")

checkpoint_files <- list.files(CONFIG$preprocessing_checkpoint, 
                               pattern = "05_final_export.*\\.rds$", 
                               full.names = TRUE)

if(length(checkpoint_files) == 0) {
  stop("ERROR: No preprocessing checkpoint found!\n",
       "Please run 01_preprocessing.R first.\n")
}

# Use most recent checkpoint
latest_checkpoint <- checkpoint_files[order(file.mtime(checkpoint_files), 
                                           decreasing = TRUE)[1]]
cat("  Using:", basename(latest_checkpoint), "\n")

preprocessed <- readRDS(latest_checkpoint)
seurat_data <- preprocessed$seurat_data
seurat_metadata <- preprocessed$seurat_metadata

cat("✓ Data loaded successfully\n")
cat("  Proteins:", nrow(seurat_data), "\n")
cat("  Cells:", ncol(seurat_data), "\n")

if("n_contaminants_removed" %in% names(preprocessed)) {
  cat("  Contaminants removed:", preprocessed$n_contaminants_removed, "\n")
}
cat("\n")

# Create Seurat object
rownames(seurat_metadata) <- seurat_metadata$Cell_ID
scp_seurat <- CreateSeuratObject(
  counts = seurat_data, 
  meta.data = seurat_metadata,
  project = "scProteomics_Harmony"
)
scp_seurat <- SetAssayData(scp_seurat, slot = "data", new.data = seurat_data)

# =============================================================================
# STEP 2: PCA
# =============================================================================

cat("=== STEP 2: PRINCIPAL COMPONENT ANALYSIS ===\n")

scp_seurat <- FindVariableFeatures(
  scp_seurat, 
  selection.method = "vst", 
  nfeatures = 1000, 
  verbose = FALSE
)

scp_seurat <- ScaleData(
  scp_seurat, 
  features = VariableFeatures(scp_seurat), 
  verbose = FALSE
)

scp_seurat <- RunPCA(
  scp_seurat, 
  features = VariableFeatures(scp_seurat), 
  npcs = 50, 
  verbose = FALSE
)

variance_explained <- (scp_seurat@reductions$pca@stdev^2) / 
  sum(scp_seurat@reductions$pca@stdev^2) * 100

cat("✓ PCA complete\n")
cat("  PC1 variance explained:", round(variance_explained[1], 1), "%\n")
cat("  PC2 variance explained:", round(variance_explained[2], 1), "%\n\n")

# =============================================================================
# STEP 3: HARMONY BATCH CORRECTION
# =============================================================================

cat("=== STEP 3: HARMONY BATCH CORRECTION ===\n")

scp_seurat <- RunHarmony(
  scp_seurat,
  group.by.vars = CONFIG$batch_variable,
  dims.use = 1:CONFIG$n_pcs,
  theta = CONFIG$theta,
  sigma = CONFIG$sigma,
  lambda = CONFIG$lambda,
  nclust = CONFIG$nclust,
  max_iter = CONFIG$max_iter,
  early.stop = TRUE,
  plot_convergence = FALSE,
  verbose = FALSE
)

cat("✓ Harmony integration complete\n\n")

# =============================================================================
# STEP 4: UMAP VISUALIZATION
# =============================================================================

cat("=== STEP 4: UMAP DIMENSIONALITY REDUCTION ===\n")

scp_seurat <- RunUMAP(
  scp_seurat,
  reduction = "harmony",
  dims = 1:CONFIG$n_pcs,
  n.neighbors = CONFIG$n_neighbors,
  min.dist = CONFIG$min_dist,
  spread = CONFIG$spread,
  metric = CONFIG$metric,
  verbose = FALSE
)

cat("✓ UMAP complete\n\n")

# =============================================================================
# STEP 5: CLUSTERING
# =============================================================================

cat("=== STEP 5: CLUSTERING ===\n")

scp_seurat <- FindNeighbors(
  scp_seurat, 
  reduction = "harmony", 
  dims = 1:CONFIG$n_pcs, 
  verbose = FALSE
)

scp_seurat <- FindClusters(
  scp_seurat, 
  resolution = CONFIG$target_resolution, 
  verbose = FALSE
)

# Rename clusters (1-indexed instead of 0-indexed)
cluster_col <- paste0("RNA_snn_res.", CONFIG$target_resolution)
original_clusters <- scp_seurat@meta.data[[cluster_col]]
renamed_clusters <- factor(as.numeric(as.character(original_clusters)) + 1)
scp_seurat@meta.data[[cluster_col]] <- renamed_clusters
Idents(scp_seurat) <- renamed_clusters
scp_seurat$seurat_clusters <- Idents(scp_seurat)

n_clusters <- length(unique(Idents(scp_seurat)))

cat("✓ Clustering complete\n")
cat("  Resolution:", CONFIG$target_resolution, "\n")
cat("  Number of clusters:", n_clusters, "\n")
cat("  Cluster IDs:", paste(sort(unique(as.character(Idents(scp_seurat)))), 
                           collapse = ", "), "\n\n")

# =============================================================================
# STEP 6: CLUSTER ENRICHMENT ANALYSIS
# =============================================================================

cat("=== STEP 6: CLUSTER ENRICHMENT ANALYSIS ===\n")

enrichment_res <- data.frame()

for(cl in unique(Idents(scp_seurat))) {
  for(grp in unique(scp_seurat$Group)) {
    in_cluster_in_group <- sum(Idents(scp_seurat) == cl & scp_seurat$Group == grp)
    in_cluster_not_group <- sum(Idents(scp_seurat) == cl & scp_seurat$Group != grp)
    not_cluster_in_group <- sum(Idents(scp_seurat) != cl & scp_seurat$Group == grp)
    not_cluster_not_group <- sum(Idents(scp_seurat) != cl & scp_seurat$Group != grp)
    
    contingency_2x2 <- matrix(
      c(in_cluster_in_group, in_cluster_not_group,
        not_cluster_in_group, not_cluster_not_group), 
      nrow = 2, byrow = TRUE
    )
    
    test_result <- fisher.test(contingency_2x2)
    
    total_in_cluster <- sum(Idents(scp_seurat) == cl)
    prop_in_cluster <- in_cluster_in_group / total_in_cluster
    total_in_group <- sum(scp_seurat$Group == grp)
    prop_in_group <- in_cluster_in_group / total_in_group
    
    enrichment_res <- rbind(enrichment_res, data.frame(
      cluster = cl,
      group = grp,
      n_cells = in_cluster_in_group,
      prop_cluster = prop_in_cluster,
      prop_group = prop_in_group,
      odds_ratio = test_result$estimate,
      p_value = test_result$p.value,
      fold_enrichment = (in_cluster_in_group / total_in_cluster) / 
        (total_in_group / ncol(scp_seurat))
    ))
  }
}

enrichment_res$p_adj <- p.adjust(enrichment_res$p_value, method = "BH")
sig_enrichments <- enrichment_res[enrichment_res$p_adj < 0.05, ]

cat("✓ Enrichment analysis complete\n")
cat("  Total tests:", nrow(enrichment_res), "\n")
cat("  Significant enrichments (FDR < 0.05):", nrow(sig_enrichments), "\n\n")

# =============================================================================
# STEP 7: VISUALIZATION
# =============================================================================

cat("=== STEP 7: CREATING VISUALIZATIONS ===\n")

# Theme for publication-quality plots
pub_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22)
  )

# Color palettes
cluster_colors <- setNames(
  colorRampPalette(c("#fee5d9", "#fc9272", "#de2d26", "#a50f15"))(n_clusters),
  sort(unique(Idents(scp_seurat)))
)

group_colors <- setNames(
  pal(scp_seurat$Group),
  sort(unique(scp_seurat$Group))
)

# Plot 1: Clusters
p_clusters <- DimPlot(
  scp_seurat, 
  reduction = "umap", 
  cols = cluster_colors, 
  pt.size = 3
) +
  ggtitle("Cluster Assignment") + 
  pub_theme

# Plot 2: Groups
p_groups <- DimPlot(
  scp_seurat, 
  reduction = "umap", 
  group.by = "Group",
  cols = group_colors, 
  pt.size = 3
) +
  ggtitle("Experimental Groups") + 
  pub_theme

# Plot 3: Batches
n_batches <- length(unique(scp_seurat@meta.data$Line_Treatment))
batch_colors <- setNames(
  hue_pal()(n_batches),
  sort(unique(scp_seurat@meta.data$Line_Treatment))
)

p_batches <- DimPlot(
  scp_seurat, 
  reduction = "umap", 
  group.by = "Line_Treatment",
  cols = batch_colors, 
  pt.size = 3
) +
  ggtitle("Batch Integration") + 
  pub_theme

# Plot 4: Treatment
p_treatment <- DimPlot(
  scp_seurat, 
  reduction = "umap", 
  group.by = "Treatment", 
  pt.size = 3
) +
  ggtitle("Treatment Groups") + 
  pub_theme

# Plot 5: Plate batch effect
p_batch <- DimPlot(
  scp_seurat, 
  reduction = "umap", 
  group.by = "Batch", 
  pt.size = 3
) +
  ggtitle("Plate Batch Effect") + 
  pub_theme

# Save plots
plots_dir <- file.path(CONFIG$output_dir, "plots")

ggsave(file.path(plots_dir, "umap_clusters.png"), p_clusters, 
       width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(plots_dir, "umap_groups.png"), p_groups, 
       width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(plots_dir, "umap_line_treatment.png"), p_batches, 
       width = 16, height = 12, dpi = 300, bg = "white")
ggsave(file.path(plots_dir, "umap_treatment.png"), p_treatment, 
       width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(plots_dir, "umap_batch.png"), p_batch, 
       width = 14, height = 12, dpi = 300, bg = "white")

# Save PDF versions
ggsave(file.path(plots_dir, "umap_clusters.pdf"), p_clusters, 
       width = 14, height = 12)
ggsave(file.path(plots_dir, "umap_groups.pdf"), p_groups, 
       width = 14, height = 12)
ggsave(file.path(plots_dir, "umap_line_treatment.pdf"), p_batches, 
       width = 16, height = 12)
ggsave(file.path(plots_dir, "umap_treatment.pdf"), p_treatment, 
       width = 14, height = 12)
ggsave(file.path(plots_dir, "umap_batch.pdf"), p_batch, 
       width = 14, height = 12)

cat("✓ UMAP plots saved\n")

# Group distribution plot
cluster_composition <- table(Cluster = Idents(scp_seurat), Group = scp_seurat$Group)
comp_proportions <- prop.table(cluster_composition, margin = 2)
comp_data <- as.data.frame(comp_proportions)
colnames(comp_data) <- c("Cluster", "Group", "Proportion")

p_group_distribution <- ggplot(comp_data, aes(x = Cluster, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Distribution of Each Group Across Clusters",
    subtitle = "Proportion of cells from each group in each cluster",
    x = "Cluster", 
    y = "Proportion of Cells"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 23),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(plots_dir, "group_distribution.png"),
  p_group_distribution, 
  width = 14, 
  height = 9, 
  dpi = 300, 
  bg = "white"
)

cat("✓ Group distribution plot saved\n\n")

# =============================================================================
# STEP 8: SAVE RESULTS
# =============================================================================

cat("=== STEP 8: SAVING RESULTS ===\n")

# Save Seurat object
saveRDS(scp_seurat, file.path(CONFIG$output_dir, "seurat_harmony.rds"))
cat("✓ Seurat object saved\n")

# Save cluster assignments
cluster_assignments <- data.frame(
  Cell_ID = colnames(scp_seurat),
  Cluster = as.character(Idents(scp_seurat)),
  Group = scp_seurat$Group,
  Line_Treatment = scp_seurat$Line_Treatment,
  stringsAsFactors = FALSE
)
write.csv(
  cluster_assignments,
  file.path(CONFIG$output_dir, "tables", "cluster_assignments.csv"),
  row.names = FALSE
)
cat("✓ Cluster assignments saved\n")

# Save enrichment results
write.csv(
  enrichment_res,
  file.path(CONFIG$output_dir, "tables", "enrichment_analysis.csv"),
  row.names = FALSE
)
cat("✓ Enrichment analysis saved\n")

# Save group distribution table
comp_table <- as.data.frame.matrix(comp_proportions)
comp_table$Cluster <- rownames(comp_table)
comp_table <- comp_table[, c("Cluster", colnames(comp_proportions))]
write.csv(
  comp_table,
  file.path(CONFIG$output_dir, "tables", "group_distribution_by_cluster.csv"),
  row.names = FALSE
)
cat("✓ Group distribution table saved\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("HARMONY INTEGRATION COMPLETE\n")
cat("=============================================================================\n\n")
cat("Summary:\n")
cat("  Proteins:", nrow(scp_seurat), "\n")
cat("  Cells:", ncol(scp_seurat), "\n")
cat("  Clusters:", n_clusters, "\n")
cat("  Significant enrichments:", nrow(sig_enrichments), "\n\n")
cat("Output directory:", CONFIG$output_dir, "\n")
cat("  ├── seurat_harmony.rds\n")
cat("  ├── plots/\n")
cat("  │   ├── umap_clusters.png\n")
cat("  │   ├── umap_groups.png\n")
cat("  │   ├── umap_line_treatment.png\n")
cat("  │   ├── umap_treatment.png\n")
cat("  │   ├── umap_batch.png\n")
cat("  │   └── group_distribution.png\n")
cat("  └── tables/\n")
cat("      ├── cluster_assignments.csv\n")
cat("      ├── enrichment_analysis.csv\n")
cat("      └── group_distribution_by_cluster.csv\n\n")
cat("Next step: Run 03_differential_analysis.R\n\n")
