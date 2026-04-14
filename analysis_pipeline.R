############################################################
# Single-cell RNA-seq analysis pipeline (MC38 tumor model)

This script contains the main steps used in the analysis:
 - Data loading and preprocessing
 - Quality control
 - Metadata integration
 - Pseudobulk differential expression (DESeq2)
 - Visualization and summary statistics

# Data are available in GEO: GSE325953, GSE326061
############################################################
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(DESeq2)
})

############################################################
# 1. Load Seurat object
# Raw data were processed using Cell Ranger (10x Genomics)
# Count matrices were imported into Seurat objects
############################################################

seurat_object <- readRDS("path_to_seurat_object.rds")

############################################################
# 2. Define experimental design
############################################################

# Treatment groups (4 groups, 12 mice total)
Bio_vec <- c(
  "Untreated",
  "N297A (Agonist, Fc-null)",
  "GA-aFuc (Agonist, Fc enhanced)",
  "GITR/Syn GA-aFuc (Non-agonist, Fc enhanced)"
)

# Mouse IDs per group
mouse_ids <- c(
  "M2","M6","M10",   # Untreated
  "M1","M4","M12",   # N297A
  "M3","M7","M11",   # GA-aFuc
  "M5","M8","M9"     # GITR/Syn
)

# Cell types
CT_vec <- c(
  "Monocytes","Macrophages","cDC1","cDC2",
  "pDC","mregDC","NK","CD8","CD4","Tregs"
)

############################################################
# 3. QC metrics
############################################################

seurat_object[["percent.mt"]] <- PercentageFeatureSet(
  seurat_object, pattern = "^mt-"
)

VlnPlot(seurat_object,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3)

############################################################
# 4. Basic visualization
############################################################

DimPlot(seurat_object, reduction = "umap", label = TRUE)

############################################################
# 5. Proportion analysis (cell composition)
############################################################

df_prop <- seurat_object@meta.data %>%
  group_by(Bio, CT) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Bio) %>%
  mutate(percent = 100 * n_cells / sum(n_cells))

ggplot(df_prop, aes(x = CT, y = percent, fill = Bio)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################
# 6. Pseudobulk DESeq2 (example)
############################################################

# extract counts
counts <- GetAssayData(seurat_object, layer = "counts")
meta <- seurat_object@meta.data

# aggregate per mouse
pseudo_bulk <- sapply(unique(meta$No), function(m) {
  cells <- rownames(meta)[meta$No == m]
  Matrix::rowSums(counts[, cells, drop = FALSE])
})

# metadata
coldata <- unique(meta[, c("No","Bio")])
rownames(coldata) <- coldata$No

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudo_bulk),
  colData = coldata,
  design = ~ Bio
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res <- results(dds)

############################################################
# 7. Save outputs
############################################################

write.csv(as.data.frame(res), "DESeq2_results.csv")

