# --- 1. SETUP: LOAD LIBRARIES ---
library(dplyr)
library(readr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
# Packages for re-generating expression data
library(tximport)
library(edgeR)
library(limma)

# --- 2. LOAD DATA ---
# Reading the .tsv file with read_tsv()
intersect_results <- readr::read_tsv("Intersection_limma_sleuth_DE_genes.tsv")
# Load the sample metadata (this was originally a CSV)
metadata_path <- "/mnt/lustre/RDS-live/bioinformatics/proj/bsp/bsp_273/nicos/bulk_de/sample_condition_path.csv"
s2c <- readr::read_csv(metadata_path)
colnames(s2c)[1] = "sample"

# --- 3. RE-GENERATE NORMALISED EXPRESSION DATA ---
files <- file.path(s2c$path, "abundance.h5")
names(files) <- s2c$sample
txi <- tximport(files, type = "kallisto", txOut = TRUE)
y <- DGEList(txi$counts)
s2c$condition <- factor(s2c$condition, levels = c("mock", "infc"))
design <- model.matrix(~condition, data = s2c)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
v <- voom(y, design, plot = FALSE) # Voom object containing normalised log-CPM values

# --- 4. PREPARE DATA FOR HEATMAP ---
# Get the list of representative transcript IDs for our intersecting genes
top_transcript_ids <- intersect_results$target_id

# Subset the voom expression matrix to include only these top transcripts
heatmap_matrix <- v$E[rownames(v$E) %in% top_transcript_ids, ]

# Create an annotation data frame for the columns (samples)
sample_annotation <- data.frame(condition = s2c$condition)
rownames(sample_annotation) <- s2c$sample

# Define colours for the conditions
annotation_colours <- list(
  condition = c(mock = "royalblue", infc = "firebrick")
)

# --- 5. GENERATE AND SAVE HEATMAP ---
print("### Generating and saving heatmap...")
# Changed height to 7 and set show_rownames to FALSE
pdf("Heatmap_Top_DE_Genes.pdf", width = 10, height = 7)
pheatmap(heatmap_matrix,
         scale = "row", # Scale genes to Z-scores to show patterns
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = sample_annotation,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         cutree_cols = 2,
         show_rownames = FALSE, # Remove gene names from the plot
         fontsize_col = 10
)
dev.off()