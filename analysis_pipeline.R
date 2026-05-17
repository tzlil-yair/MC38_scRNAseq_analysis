############################################################
# Single-cell RNA-seq analysis pipeline (MC38 tumor model)

# This script summarizes the main downstream analyses used for single-cell RNA-seq analysis of MC38 tumor-infiltrating immune cells.

The workflow includes:
 - Loading a processed and annotated Seurat object
 - Quality-control visualization
 - SCTransform normalization workflow (general template)
 - Metadata setup and experimental design
 - UMAP visualization
 - Cell composition analysis
 - Pseudobulk differential expression analysis (DESeq2)
 - Source-data export for DotPlot-style figures
 - Pathway-enrichment visualization
 - Gene-signature scoring using UCell
 - Ligand–receptor analysis using MultiNicheNet

# Raw sequencing data were processed using Cell Ranger before downstream analysis.

Data availability:
# GEO accession numbers: GSE325953, GSE326061
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(DESeq2)
  library(tidyr)
  library(readr)
})

############################################################
# 1. Load processed Seurat object
############################################################

# Replace with local path
seurat_rds <- "path_to_processed_annotated_seurat_object.rds"

# Output directory
output_dir <- "analysis_outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load object
seurat_object <- readRDS(seurat_rds)

DefaultAssay(seurat_object) <- "RNA"

############################################################
# 2. Experimental design
############################################################

# Treatment groups
condition_levels <- c(
  "Untreated",
  "N297A (Agonist, Fc-null)",
  "GA-aFuc (Agonist, Fc enhanced)",
  "GITR/Syn GA-aFuc (Non-agonist, Fc enhanced)"
)

# Mouse IDs
mouse_levels <- c(
  "M2", "M6", "M10",    # Untreated
  "M1", "M4", "M12",    # N297A
  "M3", "M7", "M11",    # GA-aFuc
  "M5", "M8", "M9"      # GITR/Syn
)

# Broad cell-type annotations
celltype_levels <- c(
  "Monocytes",
  "Macrophages",
  "cDC1",
  "cDC2",
  "pDC",
  "mregDC",
  "NK",
  "CD8",
  "CD4",
  "Tregs"
)

# Metadata formatting
seurat_object$Bio <- factor(
  seurat_object$Bio,
  levels = condition_levels
)

seurat_object$No <- factor(
  seurat_object$No,
  levels = mouse_levels
)

seurat_object$CT <- factor(
  seurat_object$CT,
  levels = celltype_levels
)

############################################################
# 3. Quality control
############################################################

# Calculate mitochondrial percentage if absent
if (!"percent.mt" %in% colnames(seurat_object@meta.data)) {

  seurat_object[["percent.mt"]] <-
    PercentageFeatureSet(
      seurat_object,
      pattern = "^mt-"
    )
}

# QC thresholds used in the study:
# genes/cell < 7,500
# UMI counts/cell < 40,000
# mitochondrial genes < 5%

qc_plot <- VlnPlot(
  seurat_object,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "No",
  ncol = 3
)

ggsave(
  file.path(output_dir, "QC_violin_plot.pdf"),
  qc_plot,
  width = 12,
  height = 5
)

############################################################
# 4. SCTransform normalization workflow
############################################################

# The preprocessing workflow used in the study included:
#
# seurat_object <- SCTransform(
#   seurat_object,
#   vst.flavor = "v2",
#   method = "glmGamPoi",
#   vars.to.regress = "percent.mt"
# )
#
# seurat_object <- RunPCA(seurat_object)
# seurat_object <- RunUMAP(seurat_object, dims = 1:30)
# seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
# seurat_object <- FindClusters(seurat_object)
#
# Harmony integration was applied when batch correction
# was required.

############################################################
# 5. Cell-type annotation
############################################################

# Cell-type annotations were assigned based on:
# - cluster marker genes
# - canonical immune-cell markers
# - ImmGen/MSigDB immune signatures

table(seurat_object$CT)

############################################################
# 6. UMAP visualization
############################################################

if ("umap" %in% names(seurat_object@reductions)) {

  p_umap_ct <- DimPlot(
    seurat_object,
    reduction = "umap",
    group.by = "CT",
    label = TRUE
  )

  p_umap_bio <- DimPlot(
    seurat_object,
    reduction = "umap",
    group.by = "Bio"
  )

  ggsave(
    file.path(output_dir, "UMAP_by_cell_type.pdf"),
    p_umap_ct,
    width = 8,
    height = 6
  )

  ggsave(
    file.path(output_dir, "UMAP_by_treatment.pdf"),
    p_umap_bio,
    width = 8,
    height = 6
  )
}

############################################################
# 7. Cell composition analysis
############################################################

cell_composition <- seurat_object@meta.data %>%
  group_by(Bio, CT) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  group_by(Bio) %>%
  mutate(
    percent = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

write.csv(
  cell_composition,
  file.path(
    output_dir,
    "cell_composition_by_treatment.csv"
  ),
  row.names = FALSE
)

p_composition <- ggplot(
  cell_composition,
  aes(
    x = CT,
    y = percent,
    fill = Bio
  )
) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    color = "black"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = "Cell type",
    y = "Percentage of cells",
    fill = "Treatment"
  )

ggsave(
  file.path(
    output_dir,
    "cell_composition_by_treatment.pdf"
  ),
  p_composition,
  width = 9,
  height = 5
)

############################################################
# 8. Pseudobulk DESeq2 analysis
############################################################

# Differential expression analysis was performed
# separately for cDC1, cDC2 and mregDC populations.

run_pseudobulk_deseq <- function(
    seurat_object,
    cell_type,
    control_condition,
    treatment_condition,
    output_dir,
    min_total_counts = 10
) {

  message(
    "Running DESeq2 for ",
    cell_type
  )

  object_subset <- subset(
    seurat_object,
    subset =
      CT == cell_type &
      Bio %in% c(
        control_condition,
        treatment_condition
      )
  )

  object_subset$Bio <- factor(
    object_subset$Bio,
    levels = c(
      control_condition,
      treatment_condition
    )
  )

  counts_matrix <- GetAssayData(
    object_subset,
    assay = "RNA",
    layer = "counts"
  )

  metadata <- object_subset@meta.data

  selected_mice <- unique(
    as.character(metadata$No)
  )

  pseudobulk_counts <- sapply(
    selected_mice,
    function(mouse_id) {

      selected_cells <- rownames(metadata)[
        metadata$No == mouse_id
      ]

      Matrix::rowSums(
        counts_matrix[
          ,
          selected_cells,
          drop = FALSE
        ]
      )
    }
  )

  pseudobulk_counts <- as.matrix(
    pseudobulk_counts
  )

  coldata <- metadata %>%
    select(No, Bio) %>%
    distinct() %>%
    filter(
      No %in% colnames(
        pseudobulk_counts
      )
    ) %>%
    as.data.frame()

  rownames(coldata) <- as.character(
    coldata$No
  )

  coldata <- coldata[
    colnames(pseudobulk_counts),
    ,
    drop = FALSE
  ]

  dds <- DESeqDataSetFromMatrix(
    countData = round(
      pseudobulk_counts
    ),
    colData = coldata,
    design = ~ Bio
  )

  dds <- dds[
    rowSums(counts(dds)) >=
      min_total_counts,
  ]

  dds$Bio <- relevel(
    dds$Bio,
    ref = control_condition
  )

  dds <- DESeq(dds)

  res <- results(dds)

  res_df <- as.data.frame(res)

  res_df$gene <- rownames(res_df)

  safe_cell_type <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    cell_type
  )

  out_folder <- file.path(
    output_dir,
    paste0(
      "DESeq2_",
      safe_cell_type
    )
  )

  dir.create(
    out_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    res_df,
    file.path(
      out_folder,
      paste0(
        safe_cell_type,
        "_DESeq2_results.csv"
      )
    ),
    row.names = FALSE
  )

  ############################################################
  # Pathway-enrichment preparation
  ############################################################

  # Top 300 genes ranked by DESeq2 Wald statistic
  # were used for pathway enrichment analysis
  # in Metascape.

  top_up <- res_df %>%
    filter(
      !is.na(stat),
      log2FoldChange > 0
    ) %>%
    arrange(desc(stat)) %>%
    slice_head(n = 300)

  write.table(
    top_up$gene,
    file.path(
      out_folder,
      paste0(
        safe_cell_type,
        "_top300_upregulated_genes.txt"
      )
    ),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  invisible(res_df)
}

# Main comparisons used in the study
for (cell_type in c(
  "cDC1",
  "cDC2",
  "mregDC"
)) {

  run_pseudobulk_deseq(
    seurat_object = seurat_object,
    cell_type = cell_type,
    control_condition =
      "N297A (Agonist, Fc-null)",
    treatment_condition =
      "GA-aFuc (Agonist, Fc enhanced)",
    output_dir = output_dir
  )
}

############################################################
# 9. Source-data export for DotPlot figures
############################################################

export_dotplot_source_data <- function(
    seurat_object,
    gene,
    cell_types,
    conditions,
    output_file
) {

  object_subset <- subset(
    seurat_object,
    subset =
      CT %in% cell_types &
      Bio %in% conditions
  )

  object_subset$CT_Bio <- paste(
    object_subset$CT,
    object_subset$Bio,
    sep = "|"
  )

  dotplot_data <- DotPlot(
    object_subset,
    features = gene,
    group.by = "CT_Bio",
    assay = "RNA"
  )$data

  dotplot_data <- dotplot_data %>%
    separate(
      id,
      into = c(
        "Cell_Type",
        "Condition"
      ),
      sep = "\\|",
      remove = FALSE
    ) %>%
    transmute(
      Cell_Type,
      Condition,
      Gene = features.plot,
      Percent_expressing = pct.exp,
      Average_expression = avg.exp,
      Average_expression_scaled =
        avg.exp.scaled
    )

  write.csv(
    dotplot_data,
    output_file,
    row.names = FALSE
  )

  invisible(dotplot_data)
}

############################################################
# 10. Pathway enrichment visualization
############################################################

# Pathway-enrichment results were visualized
# as heatmap-style annotations and ranked
# using -log10(P-values) from Metascape.

plot_metascape_barplot <- function(
    enrichment_file,
    output_file,
    plot_title = "Selected pathways"
) {

  enrichment_table <- read_csv(
    enrichment_file,
    show_col_types = FALSE
  )

  if ("LogP" %in% colnames(enrichment_table)) {

    plot_df <- enrichment_table %>%
      mutate(
        neglog10P =
          abs(as.numeric(LogP))
      )

  } else {

    stop(
      "Metascape table must contain LogP column."
    )
  }

  p <- ggplot(
    plot_df,
    aes(
      x = neglog10P,
      y = reorder(
        Description,
        neglog10P
      ),
      fill = neglog10P
    )
  ) +
    geom_col(
      color = "black"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y =
        element_blank(),
      legend.position = "none"
    ) +
    labs(
      x = expression(-log[10](P)),
      y = NULL,
      title = plot_title
    )

  ggsave(
    output_file,
    p,
    width = 6,
    height = 5
  )
}

############################################################
# 10b. Gene-signature scoring (UCell)
############################################################
############################################################
# 10b. Gene-signature scoring (UCell)
############################################################

# UCell was used to calculate single-cell
# gene-signature enrichment scores based on
# ranked gene-expression profiles.

# Example signatures included cytotoxicity-
# associated genes in CD4 T cells.

# library(UCell)
#
# cytotoxicity_signature <- list(
#   Cytotoxicity = c(
#     "Eomes",
#     "Gzmk",
#     "Ifng",
#     "Gzmb",
#     "Gzma",
#     "Nkg7",
#     "Prf1"
#   )
# )
#
# seurat_object <- AddModuleScore_UCell(
#   seurat_object,
#   features = cytotoxicity_signature
# )
#
# The resulting UCell scores were used for
# downstream visualization and comparison
# between treatment groups.
#
# UCell package: https://github.com/carmonalab/UCell
#
# Package publication:
# Andreatta & Carmona, Computational and
# Structural Biotechnology Journal (2021)
# https://doi.org/10.1016/j.csbj.2021.06.043

############################################################
# 11. Ligand–receptor analysis (MultiNicheNet)
############################################################

# Sender cells:
# Tregs
#
# Receiver cells:
# cDC1, cDC2, pDC, mregDC,
# Macrophages, Monocytes

suppressPackageStartupMessages({
  library(multinichenetr)
  library(nichenetr)
  library(SingleCellExperiment)
  library(tibble)
})

organism <- "mouse"

options(timeout = 120)

############################################################
# Load ligand–receptor network and ligand-target matrix
############################################################

lr_network_all <- readRDS(
  url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )
)

# Gene aliases were standardized prior to analysis.

lr_network <- lr_network_all %>%
  distinct(ligand, receptor)

ligand_target_matrix <- readRDS(
  url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  )
)

############################################################
# Prepare SingleCellExperiment object
############################################################

sce <- as.SingleCellExperiment(
  seurat_object,
  assay = "RNA"
)

############################################################
# Define analysis settings
############################################################

sample_id <- "No"
group_id <- "Bio"
celltype_id <- "CT"

senders_oi <- c("Tregs")

receivers_oi <- c(
  "cDC1",
  "cDC2",
  "pDC",
  "mregDC",
  "Macrophages",
  "Monocytes"
)

############################################################
# Run MultiNicheNet analysis
############################################################

multinichenet_output <-
  multi_nichenet_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi
  )

############################################################
# Extract top ligand–receptor interactions
############################################################

prioritized_tbl <-
  get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables,
    top_n = 50,
    rank_per_group = FALSE
  )

############################################################
# 12. Export metadata
############################################################

write.csv(
  as.data.frame(
    seurat_object@meta.data
  ),
  file.path(
    output_dir,
    "seurat_object_metadata.csv"
  ),
  row.names = TRUE
)

message(
  "Pipeline completed successfully."
)
