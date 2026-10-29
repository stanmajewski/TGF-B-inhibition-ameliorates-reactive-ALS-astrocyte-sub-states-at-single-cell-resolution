# Single-Cell Mass Spectrometry Analysis Pipeline

This repository contains a complete analysis pipeline for single-cell proteomics data, including preprocessing, batch correction with Harmony, clustering, and comprehensive pathway enrichment analysis.

## Overview

The pipeline consists of three main scripts that process single-cell mass spectrometry data from raw intensity matrices through to differential expression and pathway analysis:

1. **01_preprocessing.R** - Quality control, filtering, imputation, and normalization
2. **02_harmony_clustering.R** - Batch correction, dimensionality reduction, and clustering
3. **03_differential_analysis.R** - Marker identification and pathway enrichment

## Requirements

### R Version
- R >= 4.0.0

### Required R Packages

```r
# Data manipulation
install.packages(c("dplyr", "tidyr", "readxl"))

# Single-cell analysis
install.packages("Seurat")
BiocManager::install("harmony")

# Imputation
BiocManager::install("impute")

# Visualization
install.packages(c("ggplot2", "pheatmap", "patchwork", "scales", "colorspace", "RColorBrewer"))
BiocManager::install("ComplexHeatmap")

# Pathway enrichment
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",      # Human annotations
  "org.Mm.eg.db",      # Mouse annotations (if needed)
  "enrichplot",
  "ReactomePA"
))
install.packages("ggraph")
```

## Input Data Format

A single Excel file (.xlsx) exported from Spectronaut.

- Rows: Proteins
- First 5 columns: Protein metadata (column 1 = protein/gene ID)
- Columns 6+: Cell intensity values, one column per cell

### Cell naming

Sample metadata is read directly from the column names, which must follow:

```
<Line>_<Treatment>_cell<NN>      e.g. CTRL1_U_cell01, Mut1-1_SB_cell12
```

Spectronaut prefixes/suffixes such as `[1] ` and `.d.PG.Quantity` are stripped automatically.

| Field | Values | Meaning |
|---|---|---|
| Line | CTRL1, CTRL2, CTRL6, NCRMC2, Mut1-1, Mut2-1 | iPSC line |
| Treatment | `U`, `SB` | Untreated (vehicle) or SB431542 |

Derived metadata columns:
- `Genotype`: `VCP` for lines starting with `Mut`, otherwise `CTRL`
- `Group`: `Genotype_Treatment` (e.g. `CTRL_U`, `VCP_SB`)
- `Line_Treatment`: e.g. `CTRL1_U` (used for Harmony batch correction)
- `Batch`: the cell line (used for per-batch filtering and imputation)

## Directory Structure

```
project/
├── data/
│   └── raw/
│       └── proteomics_data.xlsx
├── results/
│   ├── preprocessing/
│   ├── harmony/
│   └── markers/
└── scripts/
    ├── 01_preprocessing.R
    ├── 02_harmony_clustering.R
    └── 03_differential_analysis.R
```

## Usage

### Step 1: Preprocessing

```bash
Rscript 01_preprocessing.R \
  --input data/raw/proteomics_data.xlsx \
  --output results/preprocessing
```

**Key parameters:**
- `--min-protein-detection-global`: Minimum global protein detection rate (default: 0.10)
- `--min-protein-detection-per-group`: Minimum per-group detection rate (default: 0.20)
- `--cv-threshold`: Maximum coefficient of variation (default: 3.0)
- `--log-base`: Base for log transformation (default: 2)

**Outputs:**
- `checkpoints/05_final_export_*.rds` - Processed data for next step
- Log-transformed, filtered, imputed, and normalized data matrix
- Quality control statistics

### Step 2: Harmony Integration & Clustering

```bash
Rscript 02_harmony_clustering.R \
  --input results/preprocessing/checkpoints \
  --output results/harmony \
  --batch-variable Line_Treatment \
  --resolution 0.3
```

**Key parameters:**
- `--batch-variable`: Metadata column for batch correction (default: Line_Treatment)
- `--theta`: Harmony diversity penalty (default: 1)
- `--n-pcs`: Number of principal components (default: 30)
- `--resolution`: Clustering resolution (default: 0.3)

**Outputs:**
- `seurat_harmony.rds` - Seurat object with Harmony-corrected embeddings
- UMAP visualizations (clusters, groups, batches, treatment)
- `tables/cluster_assignments.csv` - Cell-to-cluster mapping
- `tables/enrichment_analysis.csv` - Cluster enrichment statistics

### Step 3: Differential Expression & Pathway Analysis

```bash
Rscript 03_differential_analysis.R \
  --input results/harmony/seurat_harmony.rds \
  --output results/markers \
  --organism human
```

**Key parameters:**
- `--test-method`: Statistical test (default: t)
- `--logfc-threshold`: Minimum log2 fold change (default: 1)
- `--p-value-cutoff`: P-value cutoff for enrichment (default: 0.05)
- `--organism`: human or mouse (default: human)

**Outputs:**
- `tables/all_markers.csv` - Complete differential expression results
- `tables/markers_upregulated.csv` - Upregulated markers per cluster
- `tables/markers_downregulated.csv` - Downregulated markers per cluster
- `tables/pathways/` - GO, KEGG, and Reactome enrichment results
- `plots/pathways/global_reactome_overview.png` - Top pathways visualization
## License

This analysis pipeline is provided for research use.

## Contact

For questions or issues, please open an issue on the repository.
