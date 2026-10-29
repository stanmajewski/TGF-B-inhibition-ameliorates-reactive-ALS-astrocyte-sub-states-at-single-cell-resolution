#!/usr/bin/env Rscript

# =============================================================================
# CLUSTER MARKER ANALYSIS WITH PATHWAY ENRICHMENT
# =============================================================================
# This script performs differential expression analysis and pathway enrichment
# using GO, KEGG, and Reactome databases
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type = "character", 
              default = "results/harmony/seurat_harmony.rds",
              help = "Path to Harmony Seurat object [default: %default]"),
  make_option(c("-o", "--output"), type = "character", 
              default = "results/markers",
              help = "Output directory [default: %default]"),
  make_option("--test-method", type = "character", default = "t",
              help = "Statistical test method [default: %default]"),
  make_option("--min-pct", type = "double", default = 0.25,
              help = "Minimum fraction of cells expressing gene [default: %default]"),
  make_option("--logfc-threshold", type = "double", default = 1,
              help = "Minimum log2 fold change threshold [default: %default]"),
  make_option("--p-value-cutoff", type = "double", default = 0.05,
              help = "P-value cutoff for enrichment [default: %default]"),
  make_option("--q-value-cutoff", type = "double", default = 0.2,
              help = "Q-value cutoff for enrichment [default: %default]"),
  make_option("--organism", type = "character", default = "human",
              help = "Organism (human or mouse) [default: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nCluster marker identification and pathway enrichment analysis",
  epilogue = "Example: Rscript 03_differential_analysis.R -i results/harmony/seurat_harmony.rds"
)

opt <- parse_args(opt_parser)

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
  seurat_object = opt$input,
  output_dir = opt$output,
  
  test_method = opt$`test-method`,
  min_pct = opt$`min-pct`,
  logfc_threshold = opt$`logfc-threshold`,
  top_n_markers = 10,
  
  heatmap_scale = "row",
  
  enrichment = list(
    p_value_cutoff = opt$`p-value-cutoff`,
    q_value_cutoff = opt$`q-value-cutoff`,
    min_genes = 5,
    max_genes = 500,
    organism = opt$organism
  )
)

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("CLUSTER MARKER ANALYSIS WITH PATHWAY ENRICHMENT\n")
cat("=============================================================================\n\n")

cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(colorspace)
  library(clusterProfiler)
  library(enrichplot)
  library(ggraph)
  library(ReactomePA)
  
  # Load organism database
  if(CONFIG$enrichment$organism == "human") {
    library(org.Hs.eg.db)
    org_db <- org.Hs.eg.db
  } else if(CONFIG$enrichment$organism == "mouse") {
    library(org.Mm.eg.db)
    org_db <- org.Mm.eg.db
  } else {
    stop("Organism must be 'human' or 'mouse'")
  }
})

# Create output directories
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "plots", "pathways"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "tables", "pathways"), showWarnings = FALSE, recursive = TRUE)

cat("✓ Packages loaded\n\n")

cat("Configuration:\n")
cat("  Input Seurat object:", CONFIG$seurat_object, "\n")
cat("  Output directory:", CONFIG$output_dir, "\n")
cat("  Test method:", CONFIG$test_method, "\n")
cat("  Log2FC threshold:", CONFIG$logfc_threshold, "\n")
cat("  Organism:", CONFIG$enrichment$organism, "\n\n")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

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

# Gene ID conversion function
convert_to_entrez <- function(gene_symbols) {
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  
  gene_ids <- tryCatch({
    AnnotationDbi::select(
      org_db, 
      keys = gene_symbols,
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "SYMBOL"
    )
  }, error = function(e) {
    cat("⚠️  Error in gene ID conversion:", e$message, "\n")
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })
  
  gene_ids <- gene_ids[!is.na(gene_ids$ENTREZID), ]
  gene_ids <- gene_ids[!duplicated(gene_ids$SYMBOL), ]
  
  return(gene_ids)
}

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("=== STEP 1: LOADING SEURAT OBJECT ===\n")

scp_seurat <- readRDS(CONFIG$seurat_object)

n_clusters <- length(unique(Idents(scp_seurat)))
cluster_colors <- setNames(
  colorRampPalette(c("#fee5d9", "#fc9272", "#de2d26", "#a50f15"))(n_clusters),
  sort(unique(Idents(scp_seurat)))
)

cat("✓ Data loaded\n")
cat("  Features:", nrow(scp_seurat), "\n")
cat("  Cells:", ncol(scp_seurat), "\n")
cat("  Clusters:", n_clusters, "\n\n")

# =============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

cat("=== STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
cat("Finding both upregulated and downregulated markers...\n\n")

all_markers <- data.frame()

for(cluster_id in sort(unique(Idents(scp_seurat)))) {
  cat("  Processing Cluster", cluster_id, "...\n")
  
  markers <- FindMarkers(
    scp_seurat,
    ident.1 = cluster_id,
    ident.2 = NULL,
    test.use = CONFIG$test_method,
    min.pct = CONFIG$min_pct,
    logfc.threshold = CONFIG$logfc_threshold,
    only.pos = FALSE,
    verbose = FALSE
  )
  
  markers$cluster <- cluster_id
  markers$gene <- rownames(markers)
  markers$direction <- ifelse(markers$avg_log2FC > 0, "Up", "Down")
  
  all_markers <- rbind(all_markers, markers)
  
  n_up <- sum(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05)
  n_down <- sum(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05)
  
  cat("    ↑", n_up, "upregulated | ↓", n_down, "downregulated (FDR < 0.05)\n")
}

all_markers$abs_avg_log2FC <- abs(all_markers$avg_log2FC)

markers_up <- all_markers %>%
  filter(avg_log2FC > 0) %>%
  arrange(cluster, desc(avg_log2FC))

markers_down <- all_markers %>%
  filter(avg_log2FC < 0) %>%
  arrange(cluster, avg_log2FC)

cat("\n✓ Differential expression complete\n")
cat("  Total upregulated markers:", nrow(markers_up), "\n")
cat("  Total downregulated markers:", nrow(markers_down), "\n\n")

# Save DE results
write.csv(
  all_markers,
  file.path(CONFIG$output_dir, "tables", "all_markers.csv"),
  row.names = FALSE
)
write.csv(
  markers_up,
  file.path(CONFIG$output_dir, "tables", "markers_upregulated.csv"),
  row.names = FALSE
)
write.csv(
  markers_down,
  file.path(CONFIG$output_dir, "tables", "markers_downregulated.csv"),
  row.names = FALSE
)

# =============================================================================
# STEP 3: GENE ID CONVERSION
# =============================================================================

cat("=== STEP 3: GENE ID CONVERSION ===\n")

all_gene_ids <- convert_to_entrez(all_markers$gene)

cat("✓ Gene ID conversion complete\n")
cat("  Input genes:", length(unique(all_markers$gene)), "\n")
cat("  Converted to Entrez IDs:", nrow(all_gene_ids), "\n")
cat("  Conversion rate:", 
    round(nrow(all_gene_ids) / length(unique(all_markers$gene)) * 100, 1), "%\n\n")

# =============================================================================
# STEP 4: PATHWAY ENRICHMENT ANALYSIS
# =============================================================================

cat("=== STEP 4: PATHWAY ENRICHMENT ANALYSIS ===\n\n")

all_go_results <- list()
all_kegg_results <- list()
all_reactome_results <- list()
enrichment_summary <- data.frame()

for(cluster_id in sort(unique(Idents(scp_seurat)))) {
  cat("Analyzing Cluster", cluster_id, "...\n")
  
  # Get significant genes
  cluster_genes_up <- all_markers %>%
    filter(
      cluster == cluster_id, 
      avg_log2FC > 0,
      p_val_adj < CONFIG$enrichment$p_value_cutoff
    )
  
  cluster_genes_down <- all_markers %>%
    filter(
      cluster == cluster_id,
      avg_log2FC < 0,
      p_val_adj < CONFIG$enrichment$p_value_cutoff
    )
  
  cat("  ↑", nrow(cluster_genes_up), "genes | ↓", nrow(cluster_genes_down), "genes\n")
  
  # Convert to Entrez IDs
  genes_up_entrez <- all_gene_ids %>%
    filter(SYMBOL %in% cluster_genes_up$gene) %>%
    pull(ENTREZID)
  
  genes_down_entrez <- all_gene_ids %>%
    filter(SYMBOL %in% cluster_genes_down$gene) %>%
    pull(ENTREZID)
  
  # GO enrichment (upregulated)
  if(length(genes_up_entrez) >= CONFIG$enrichment$min_genes) {
    go_up <- tryCatch({
      enrichGO(
        gene = genes_up_entrez,
        OrgDb = org_db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = CONFIG$enrichment$p_value_cutoff,
        qvalueCutoff = CONFIG$enrichment$q_value_cutoff,
        readable = TRUE
      )
    }, error = function(e) NULL)
    
    if(!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
      all_go_results[[paste0("cluster", cluster_id, "_up")]] <- as.data.frame(go_up)
      cat("    GO (↑):", nrow(as.data.frame(go_up)), "terms\n")
    }
  }
  
  # GO enrichment (downregulated)
  if(length(genes_down_entrez) >= CONFIG$enrichment$min_genes) {
    go_down <- tryCatch({
      enrichGO(
        gene = genes_down_entrez,
        OrgDb = org_db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = CONFIG$enrichment$p_value_cutoff,
        qvalueCutoff = CONFIG$enrichment$q_value_cutoff,
        readable = TRUE
      )
    }, error = function(e) NULL)
    
    if(!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
      all_go_results[[paste0("cluster", cluster_id, "_down")]] <- as.data.frame(go_down)
      cat("    GO (↓):", nrow(as.data.frame(go_down)), "terms\n")
    }
  }
  
  # KEGG enrichment (upregulated)
  if(length(genes_up_entrez) >= CONFIG$enrichment$min_genes) {
    kegg_up <- tryCatch({
      enrichKEGG(
        gene = genes_up_entrez,
        organism = ifelse(CONFIG$enrichment$organism == "human", "hsa", "mmu"),
        pAdjustMethod = "BH",
        pvalueCutoff = CONFIG$enrichment$p_value_cutoff,
        qvalueCutoff = CONFIG$enrichment$q_value_cutoff
      )
    }, error = function(e) NULL)
    
    if(!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      all_kegg_results[[paste0("cluster", cluster_id, "_up")]] <- as.data.frame(kegg_up)
      cat("    KEGG (↑):", nrow(as.data.frame(kegg_up)), "pathways\n")
    }
  }
  
  # KEGG enrichment (downregulated)
  if(length(genes_down_entrez) >= CONFIG$enrichment$min_genes) {
    kegg_down <- tryCatch({
      enrichKEGG(
        gene = genes_down_entrez,
        organism = ifelse(CONFIG$enrichment$organism == "human", "hsa", "mmu"),
        pAdjustMethod = "BH",
        pvalueCutoff = CONFIG$enrichment$p_value_cutoff,
        qvalueCutoff = CONFIG$enrichment$q_value_cutoff
      )
    }, error = function(e) NULL)
    
    if(!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      all_kegg_results[[paste0("cluster", cluster_id, "_down")]] <- as.data.frame(kegg_down)
      cat("    KEGG (↓):", nrow(as.data.frame(kegg_down)), "pathways\n")
    }
  }
  
  # Reactome enrichment (upregulated)
  if(length(genes_up_entrez) >= CONFIG$enrichment$min_genes) {
    reactome_up <- tryCatch({
      enrichPathway(
        gene = genes_up_entrez,
        organism = CONFIG$enrichment$organism,
        pAdjustMethod = "BH",
        pvalueCutoff = CONFIG$enrichment$p_value_cutoff,
        qvalueCutoff = CONFIG$enrichment$q_value_cutoff,
        readable = TRUE
      )
    }, error = function(e) NULL)
    
    if(!is.null(reactome_up) && nrow(as.data.frame(reactome_up)) > 0) {
      all_reactome_results[[paste0("cluster", cluster_id, "_up")]] <- as.data.frame(reactome_up)
      cat("    Reactome (↑):", nrow(as.data.frame(reactome_up)), "pathways\n")
    }
  }
  
  # Reactome enrichment (downregulated)
  if(length(genes_down_entrez) >= CONFIG$enrichment$min_genes) {
    reactome_down <- tryCatch({
      enrichPathway(
        gene = genes_down_entrez,
        organism = CONFIG$enrichment$organism,
        pAdjustMethod = "BH",
        pvalueCutoff = CONFIG$enrichment$p_value_cutoff,
        qvalueCutoff = CONFIG$enrichment$q_value_cutoff,
        readable = TRUE
      )
    }, error = function(e) NULL)
    
    if(!is.null(reactome_down) && nrow(as.data.frame(reactome_down)) > 0) {
      all_reactome_results[[paste0("cluster", cluster_id, "_down")]] <- as.data.frame(reactome_down)
      cat("    Reactome (↓):", nrow(as.data.frame(reactome_down)), "pathways\n")
    }
  }
  
  cat("\n")
}

cat("✓ Pathway enrichment complete\n\n")

# =============================================================================
# STEP 5: SAVE ENRICHMENT RESULTS
# =============================================================================

cat("=== STEP 5: SAVING ENRICHMENT RESULTS ===\n")

# Save GO results
for(name in names(all_go_results)) {
  write.csv(
    all_go_results[[name]],
    file.path(CONFIG$output_dir, "tables", "pathways", paste0("GO_", name, ".csv")),
    row.names = FALSE
  )
}

# Save KEGG results
for(name in names(all_kegg_results)) {
  write.csv(
    all_kegg_results[[name]],
    file.path(CONFIG$output_dir, "tables", "pathways", paste0("KEGG_", name, ".csv")),
    row.names = FALSE
  )
}

# Save Reactome results
for(name in names(all_reactome_results)) {
  write.csv(
    all_reactome_results[[name]],
    file.path(CONFIG$output_dir, "tables", "pathways", paste0("Reactome_", name, ".csv")),
    row.names = FALSE
  )
}

cat("✓ Enrichment results saved\n")
cat("  GO enrichments:", length(all_go_results), "\n")
cat("  KEGG enrichments:", length(all_kegg_results), "\n")
cat("  Reactome enrichments:", length(all_reactome_results), "\n\n")

# =============================================================================
# STEP 6: GLOBAL REACTOME OVERVIEW
# =============================================================================

cat("=== STEP 6: GLOBAL REACTOME OVERVIEW ===\n")

global_reactome_data <- data.frame()

for(cluster_id in sort(unique(Idents(scp_seurat)))) {
  reactome_up_key <- paste0("cluster", cluster_id, "_up")
  
  if(reactome_up_key %in% names(all_reactome_results)) {
    reactome_data <- all_reactome_results[[reactome_up_key]]
    
    if(nrow(reactome_data) > 0) {
      top3 <- reactome_data %>%
        arrange(p.adjust) %>%
        slice_head(n = 3) %>%
        mutate(
          cluster = as.character(cluster_id),
          rank = 1:n(),
          neg_log_p = -log10(p.adjust)
        ) %>%
        select(cluster, rank, Description, Count, p.adjust, neg_log_p, GeneRatio)
      
      global_reactome_data <- rbind(global_reactome_data, top3)
    }
  }
}

if(nrow(global_reactome_data) > 0) {
  global_reactome_data <- global_reactome_data %>%
    mutate(
      Description_short = ifelse(
        nchar(Description) > 45, 
        paste0(substr(Description, 1, 42), "..."), 
        Description
      ),
      plot_order = as.numeric(factor(cluster)) * 10 + rank
    ) %>%
    arrange(desc(plot_order))
  
  global_reactome_data$plot_label <- factor(
    paste0("C", global_reactome_data$cluster, " (", global_reactome_data$rank, 
           "): ", global_reactome_data$Description_short),
    levels = paste0("C", global_reactome_data$cluster, " (", global_reactome_data$rank, 
                   "): ", global_reactome_data$Description_short)
  )
  
  p_global_reactome <- ggplot(
    global_reactome_data, 
    aes(x = factor(cluster), y = Description_short)
  ) +
    geom_point(
      aes(size = Count, fill = neg_log_p), 
      shape = 21, color = "black", stroke = 0.8, alpha = 0.9
    ) +
    scale_size_continuous(name = "Gene\nCount", range = c(6, 20)) +
    scale_fill_gradient2(
      low = "#f7f7f7", 
      mid = "#fdae61",
      high = "#d73027",
      midpoint = median(global_reactome_data$neg_log_p),
      name = "-log10\n(adj p-val)"
    ) +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
      title = "Top 3 Upregulated Reactome Pathways per Cluster",
      x = "Cluster",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      strip.text.y.left = element_text(angle = 0, size = 20, face = "bold", hjust = 1),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95", color = NA),
      axis.text.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  ggsave(
    file.path(CONFIG$output_dir, "plots", "pathways", "global_reactome_overview.png"),
    p_global_reactome, 
    width = 15,
    height = max(10, length(unique(global_reactome_data$cluster)) * 2),
    dpi = 300, 
    bg = "white"
  )
  
  write.csv(
    global_reactome_data,
    file.path(CONFIG$output_dir, "tables", "pathways", "global_reactome_top3.csv"),
    row.names = FALSE
  )
  
  cat("✓ Global Reactome overview created\n")
  cat("  Total pathways displayed:", nrow(global_reactome_data), "\n\n")
} else {
  cat("⚠️  No Reactome pathways available\n\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("DIFFERENTIAL ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")
cat("Summary:\n")
cat("  Clusters analyzed:", n_clusters, "\n")
cat("  Upregulated markers:", nrow(markers_up), "\n")
cat("  Downregulated markers:", nrow(markers_down), "\n")
cat("  GO enrichments:", length(all_go_results), "\n")
cat("  KEGG enrichments:", length(all_kegg_results), "\n")
cat("  Reactome enrichments:", length(all_reactome_results), "\n\n")
cat("Output directory:", CONFIG$output_dir, "\n")
cat("  ├── tables/\n")
cat("  │   ├── all_markers.csv\n")
cat("  │   ├── markers_upregulated.csv\n")
cat("  │   ├── markers_downregulated.csv\n")
cat("  │   └── pathways/\n")
cat("  │       ├── GO_cluster*_*.csv\n")
cat("  │       ├── KEGG_cluster*_*.csv\n")
cat("  │       └── Reactome_cluster*_*.csv\n")
cat("  └── plots/\n")
cat("      └── pathways/\n")
cat("          └── global_reactome_overview.png\n\n")
cat("Analysis pipeline complete!\n\n")
