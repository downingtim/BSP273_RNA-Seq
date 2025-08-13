##########
# BSP273: Differential Expression Analysis Script (Final Revised Version)
# Workflow: Kallisto -> tximport -> limma-voom & sleuth
# Includes EDA, direct BioMart query, robust gene-level summarisation, DE analysis, 
# visualisation, and corrected multi-file output.
##########


#################################################################
### 1. SETUP: LOAD LIBRARIES
#################################################################

# Note: Run these install commands once if packages are not already installed.
# if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
# BiocManager::install(c("rhdf5", "limma", "edgeR", "tximport"))
# install.packages(c("devtools", "ggplot2", "ggrepel", "dplyr", "readr", "tibble", "httr"))
# devtools::install_github("pachterlab/sleuth")

# Data wrangling and plotting
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggrepel)
library(httr) # For direct web queries

# Differential expression packages
library(tximport)
library(edgeR)
library(limma)
library(sleuth)


#################################################################
### 2. DATA LOADING AND PREPARATION
#################################################################

metadata_path <- "/mnt/lustre/RDS-live/bioinformatics/proj/bsp/bsp_273/nicos/bulk_de/sample_condition_path.csv"
s2c <- readr::read_csv(metadata_path)
colnames(s2c)[1] = "sample"
print("### Sample Metadata:")
print(s2c)

so <- sleuth_prep(s2c,
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm = TRUE)


#################################################################
### 3. EXPLORATORY DATA ANALYSIS (EDA)
#################################################################

# --- 3.1. Principal Component Analysis (PCA) ---
tpm_matrix <- so$obs_norm %>%
  dplyr::select(target_id, sample, tpm) %>%
  tidyr::pivot_wider(names_from = sample, values_from = tpm) %>%
  tibble::column_to_rownames(var = "target_id") %>%
  as.matrix()

tpm_matrix_filtered <- tpm_matrix[apply(tpm_matrix, 1, var) > 0, ]
pca_result <- prcomp(t(tpm_matrix_filtered), scale. = TRUE)

pca_df <- as.data.frame(pca_result$x) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::left_join(so$sample_to_covariates, by = "sample")

var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pdf("PCA_Infected_vs_Mock.pdf", width = 7, height = 5)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.8) +
  labs(x = paste0("PC1 (", var_explained[1], "% variance)"),
       y = paste0("PC2 (", var_explained[2], "% variance)"),
       colour = "Condition", shape = "Condition") +
  theme_bw(base_size = 14)
dev.off()


# --- 3.2. Sample Correlation Heatmap ---
pdf("Sample_Correlation_Heatmap.pdf", width = 8, height = 8)
plot_sample_heatmap(so, use_filtered = TRUE)
dev.off()


#################################################################
### 4. GENE ANNOTATION FROM BIOMART (ROBUST METHOD)
#################################################################

print("### Querying Ensembl BioMart server directly...")
host <- "https://oct2024.archive.ensembl.org"
query_xml <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name="ggallus_gene_ensembl" interface="default">
        <Attribute name="ensembl_transcript_id" />
        <Attribute name="transcript_version" />
        <Attribute name="ensembl_gene_id" />
        <Attribute name="external_gene_name" />
        <Attribute name="description" />
    </Dataset>
</Query>'

response <- httr::POST(paste0(host, "/biomart/martservice"), body = list(query = query_xml))
httr::stop_for_status(response)

t2g <- readr::read_tsv(httr::content(response, "text", encoding = "UTF-8"),
                       col_names = c("ensembl_transcript_id", "transcript_version",
                                     "ensembl_gene_id", "external_gene_name", "description"))

t2g <- t2g %>%
  dplyr::rename(ens_gene = ensembl_gene_id,
                ext_gene = external_gene_name) %>%
  dplyr::mutate(target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))

# ** FIX FOR DUPLICATES **
# Ensure every transcript ID is unique to prevent duplication during joins.
t2g <- t2g %>%
  dplyr::distinct(target_id, .keep_all = TRUE)

print("### Annotation table retrieved and de-duplicated successfully:")
head(t2g)


#################################################################
### 5. DIFFERENTIAL EXPRESSION WITH LIMMA-VOOM
#################################################################

files <- file.path(s2c$path, "abundance.h5")
names(files) <- s2c$sample
txi <- tximport(files, type = "kallisto", txOut = TRUE)
y <- DGEList(txi$counts)
s2c$condition <- factor(s2c$condition, levels = c("mock", "infc"))
design <- model.matrix(~condition, data = s2c)

keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
v <- voom(y, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)
results_limma_transcript <- topTable(fit, coef = "conditioninfc", n = Inf)

annotated_transcripts <- results_limma_transcript %>%
  tibble::rownames_to_column(var = "target_id") %>%
  dplyr::left_join(t2g, by = "target_id")

# Create a robust gene-level summary
gene_level_limma_results <- annotated_transcripts %>%
  dplyr::filter(!is.na(ens_gene)) %>%
  dplyr::group_by(ens_gene) %>%
  dplyr::slice_min(order_by = adj.P.Val, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  # ** FIX FOR DUPLICATES **: Ensure unique gene IDs
  dplyr::distinct(ens_gene, .keep_all = TRUE)

# Add significance flags and labels
gene_level_limma_results <- gene_level_limma_results %>%
  dplyr::mutate(
    gene_label = case_when(!is.na(ext_gene) & ext_gene != "" ~ ext_gene, TRUE ~ ens_gene),
    significant = case_when(
      adj.P.Val < 0.05 & logFC > 2  ~ "Upregulated in Infected",
      adj.P.Val < 0.05 & logFC < -2 ~ "Downregulated in Infected",
      TRUE                         ~ "Not Significant"
    )
  )

# --- 5.5. Save Limma Output Files ---
print("### Saving gene-level limma results as TSV files...")
# File 1: All genes from limma analysis
readr::write_tsv(gene_level_limma_results, "limma_gene_level_all.tsv")

# File 2: DE genes only from limma analysis
limma_de_genes_only <- gene_level_limma_results %>%
  dplyr::filter(significant != "Not Significant")
readr::write_tsv(limma_de_genes_only, "limma_gene_level_DE_only.tsv")


# --- 5.6. Visualise Gene-Level Limma Results ---
pdf("Volcano_Plot_limma_gene_level.pdf", width = 8, height = 6)
ggplot(gene_level_limma_results, aes(x = logFC, y = -log10(adj.P.Val), colour = significant)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_colour_manual(name = "Significance",
                      values = c("Upregulated in Infected" = "red", "Downregulated in Infected" = "blue", "Not Significant" = "grey")) +
  geom_text_repel(data = head(limma_de_genes_only, 15),
                  aes(label = gene_label), size = 3.5, box.padding = 0.5, max.overlaps = Inf) +
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)") +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")
dev.off()


#################################################################
### 6. GENE-LEVEL DE ANALYSIS WITH SLEUTH
#################################################################

so_gene <- sleuth_prep(s2c, target_mapping = t2g, aggregation_column = 'ens_gene',
                       extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

so_gene <- sleuth_fit(so_gene, ~condition, 'full')
so_gene <- sleuth_fit(so_gene, ~1, 'reduced')
so_gene <- sleuth_lrt(so_gene, 'reduced', 'full')
results_sleuth <- sleuth_results(so_gene, 'reduced:full', 'lrt', show_all = TRUE) %>%
  # ** FIX FOR DUPLICATES **: Ensure unique gene IDs
  dplyr::distinct(target_id, .keep_all = TRUE)

# --- 6.1. Save Sleuth Output Files ---
print("### Saving gene-level sleuth results as TSV files...")
# File 3: All genes from sleuth analysis
readr::write_tsv(results_sleuth, "sleuth_gene_level_all.tsv")

# File 4: DE genes only from sleuth analysis
sleuth_de_genes_only <- results_sleuth %>%
  dplyr::filter(qval <= 0.05)
readr::write_tsv(sleuth_de_genes_only, "sleuth_gene_level_DE_only.tsv")


#################################################################
### 7. INTERSECTION OF LIMMA AND SLEUTH RESULTS
#################################################################

# --- 7.1. Get lists of significant gene IDs from both methods ---
limma_significant_genes <- limma_de_genes_only %>% dplyr::pull(ens_gene)
sleuth_significant_genes <- sleuth_de_genes_only %>% dplyr::pull(target_id)

intersecting_genes <- intersect(limma_significant_genes, sleuth_significant_genes)
print(paste("Number of DE genes found in both Limma and Sleuth:", length(intersecting_genes)))

# --- 7.2. Prepare and join data for the intersecting genes ---
limma_intersect_data <- gene_level_limma_results %>%
  dplyr::filter(ens_gene %in% intersecting_genes)

sleuth_intersect_data <- results_sleuth %>%
  dplyr::filter(target_id %in% intersecting_genes) %>%
  dplyr::select(ens_gene = target_id, sleuth_qval = qval)

combined_intersect_df <- dplyr::inner_join(limma_intersect_data, sleuth_intersect_data, by = "ens_gene")

# --- 7.3. Save the intersecting dataset to a TSV file ---
# File 5: The high-confidence intersection of DE genes
readr::write_tsv(combined_intersect_df, "Intersection_limma_sleuth_DE_genes.tsv")


# --- 7.4. Create a Volcano Plot for the intersecting genes only ---
pdf("Volcano_Plot_Intersection.pdf", width = 8, height = 6)
ggplot(combined_intersect_df, aes(x = logFC, y = -log10(adj.P.Val), colour = significant)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(name = "Direction",
                      values = c("Upregulated in Infected" = "red", "Downregulated in Infected" = "blue")) +
  geom_text_repel(aes(label = gene_label), size = 3.5, box.padding = 0.5, max.overlaps = 20) +
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)") +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")
dev.off()

print("### Analysis Complete. All files saved.")