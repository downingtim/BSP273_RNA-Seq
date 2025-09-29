# BSP273_RNA-Seq

Differential gene expression analysis of BSP273 data 

# Differential Expression Analysis of Infected vs. Mock Samples

This repository contains the scripts and results for a differential 
expression (DE) analysis of RNA-Seq data, comparing infected samples 
against mock controls. The analysis identifies significantly 
differentially expressed genes using two robust statistical methods 
(limma-voom and sleuth) and visualises the results.


Setup and Dependencies
======================

This analysis was performed using R v4.2.2 and requires a specific 
Conda environment to ensure all package dependencies are compatible.


1. Create Conda Environment
First, create and activate the Conda environment. This only needs to 
be done once.

2. Create the environment with all necessary R packages

    conda create -n r_de_env -c conda-forge -c bioconda \
    r-base=4.2.2 r-essentials bioconductor-limma \
    bioconductor-edger bioconductor-tximport \
    bioconductor-sleuth r-pheatmap r-httr r-ggplot2 r-ggrepel

   Activate the environment:   conda activate r_de_env


3. Required Data

Place the following file in the root of the project directory:
  * sample_condition_path.csv: A CSV file containing the sample 
    metadata with three columns: sample, condition, and path.


4. Analysis Workflow

The analysis is performed in two main steps using the provided R 
scripts. Ensure you have activated the Conda environment 
(conda activate r_de_env) before running them.

4.1: Differential Expression Analysis

The DE.R script is the primary analysis pipeline. It performs the 
end-to-end differential expression analysis.

To run:

  Rscript DE.R

This script will:

  * Load sample metadata and Kallisto quantification results.
  * Perform exploratory data analysis, generating a PCA plot.
  * Query the Ensembl BioMart database (Oct 2024 archive) to 
    annotate transcripts with gene names.
  * Run a transcript-level DE analysis using the limma-voom pipeline.
  * Run a gene-level DE analysis using sleuth.
  * Summarise the limma results at the gene level.
  * Identify the intersection of significant genes found by both 
    methods.
  * Generate all primary output files (see below).


4.2: Heatmap Visualisation

The heat.R script takes the high-confidence intersection results 
from the first step and creates a clustered heatmap.

To run:

  Rscript heat.R

This script will:

  * Load the intersection results file 
    (Intersection_limma_sleuth_DE_genes.tsv).
  * Re-generate the normalised expression data from the limma-voom 
    analysis.
  * Subset the expression matrix to include only the 375 
    high-confidence genes.
  * Generate a clustered heatmap, scaling by gene (row) to visualise 
    expression patterns.


4.3 Output Files

Running these scripts will generate the following output files in 
the project directory:


Plots (.pdf)
------------

  * PCA_Infected_vs_Mock.pdf
    PCA plot showing sample clustering.

  * Volcano_Plot_limma_gene_level.pdf
    Volcano plot of all genes from the limma analysis.

  * Volcano_Plot_Intersection.pdf
    Volcano plot showing only the high-confidence intersecting genes.

  * Heatmap_Top_DE_Genes.pdf
    Clustered heatmap of the intersecting genes.


Data Tables (.tsv)
------------------

  * limma_gene_level_all.tsv
    Gene-level results for all genes from the limma analysis.

  * limma_gene_level_DE_only.tsv
    A filtered table of only the significant DE genes from limma.

  * sleuth_gene_level_all.tsv
    Gene-level results for all genes from the sleuth analysis.

  * sleuth_gene_level_DE_only.tsv
    A filtered table of only the significant DE genes from sleuth.

  * Intersection_limma_sleuth_DE_genes.tsv
    The final table of 333 high-confidence genes found by both 
    methods, with their associated statistics.
