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

A single Excel file (.xlsx) exported from Spectronaut. No separate metadata file is needed.

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

## Analysis Methods

### Preprocessing Pipeline

1. **Data Loading**: Loads raw intensity data and parses sample metadata from cell names
2. **Log Transformation**: log2(x + 0.01) transformation
3. **Contaminant Removal**: Filters out common contaminants (keratin, trypsin, etc.)
4. **Quality Filtering**:
   - Global protein detection rate
   - Per-group detection rate
   - Per-batch detection rate
   - Missing value threshold
   - Coefficient of variation filter
5. **Hybrid Imputation**:
   - Low missing (<20%): KNN imputation
   - Medium missing (20-50%): Left-censored imputation
   - High missing (>50%): 5th percentile imputation
6. **Median-Center Normalization**: Centers each cell to global median

### Harmony Integration

1. **PCA**: Principal component analysis on variable features
2. **Harmony**: Batch effect correction using specified batch variable
3. **UMAP**: Dimensionality reduction for visualization
4. **Clustering**: Graph-based clustering with specified resolution
5. **Enrichment Testing**: Fisher's exact test for group enrichment in clusters

### Pathway Enrichment

1. **Differential Expression**: FindMarkers for each cluster vs all others
2. **Gene ID Conversion**: Gene symbols to Entrez IDs
3. **GO Enrichment**: Biological Process terms
4. **KEGG Enrichment**: KEGG pathway database
5. **Reactome Enrichment**: Reactome pathway database
6. **Bidirectional Analysis**: Both upregulated and downregulated genes

## Output Files

### Preprocessing
- `checkpoints/` - Intermediate data checkpoints
- Quality control statistics in stdout

### Harmony Clustering
- `seurat_harmony.rds` - Main Seurat object
- `plots/*.png` - UMAP visualizations
- `plots/*.pdf` - Publication-quality PDFs
- `tables/cluster_assignments.csv` - Cell metadata with clusters
- `tables/enrichment_analysis.csv` - Statistical enrichment tests
- `tables/group_distribution_by_cluster.csv` - Composition matrix

### Differential Analysis
- `tables/all_markers.csv` - All DE results
- `tables/markers_upregulated.csv` - Positive markers
- `tables/markers_downregulated.csv` - Negative markers
- `tables/pathways/GO_cluster*_*.csv` - GO enrichment per cluster
- `tables/pathways/KEGG_cluster*_*.csv` - KEGG enrichment per cluster
- `tables/pathways/Reactome_cluster*_*.csv` - Reactome enrichment per cluster
- `plots/pathways/global_reactome_overview.png` - Global pathway overview

## Statistical Methods

- **Differential Expression**: Welch's t-test (default) or Wilcoxon rank-sum test
- **Multiple Testing Correction**: Benjamini-Hochberg FDR
- **Enrichment Testing**: Fisher's exact test for cluster composition
- **Pathway Enrichment**: Hypergeometric test with FDR correction

## Customization

All scripts support command-line arguments for parameter customization. Use `--help` flag for full options:

```bash
Rscript 01_preprocessing.R --help
Rscript 02_harmony_clustering.R --help
Rscript 03_differential_analysis.R --help
```

## Citation

If you use this pipeline in your research, please cite:

- Seurat: Hao et al. (2021) Cell
- Harmony: Korsunsky et al. (2019) Nature Methods
- clusterProfiler: Yu et al. (2012) OMICS
- ReactomePA: Yu & He (2016) Molecular BioSystems

## License

This analysis pipeline is provided for research use.

## Contact

For questions or issues, please open an issue on the repository.
