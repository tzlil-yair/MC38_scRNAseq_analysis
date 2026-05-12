# MC38 scRNA-seq analysis

This repository contains code used for the analysis of single-cell RNA sequencing (scRNA-seq) data from MC38 tumor-infiltrating immune cells.

## Data availability

The data generated in this study are available in the Gene Expression Omnibus (GEO):

- GSE325953  
- GSE326061  

## Analysis overview

The analysis was performed using a processed and annotated Seurat object.

The pipeline includes:

- Quality-control (QC) visualization
- SCTransform normalization workflow
- UMAP visualization
- Cell-type annotation
- Cell composition analysis
- Pseudobulk differential expression analysis (DESeq2)
- Pathway enrichment analysis (Metascape)
- Gene-signature scoring (UCell)
- Ligand–receptor interaction analysis (MultiNicheNet)
- Source-data export for figure generation

## Main script

analysis_pipeline.R

## Notes

Raw sequencing data were processed using Cell Ranger (10x Genomics) prior to downstream analysis.

This repository provides a simplified and reproducible version of the main analysis workflow used in the study.
