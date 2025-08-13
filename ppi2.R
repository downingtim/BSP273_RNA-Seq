# ppi_network_analysis.R (Revised)

# --- 1. SETUP: LOAD LIBRARIES ---
# BiocManager::install("STRINGdb")
library(STRINGdb)
library(dplyr)
library(readr)
library(igraph)

# --- 2. LOAD DATA AND PREPARE GENE LIST ---
print("### Loading intersection results...")
intersect_results <- readr::read_tsv("Intersection_limma_sleuth_DE_genes.tsv")

# Create a clean data frame with unique gene labels and their corresponding logFC
de_genes <- intersect_results %>%
  dplyr::select(gene_label, logFC) %>%
  dplyr::distinct(gene_label, .keep_all = TRUE)

print(paste("### Preparing to map", nrow(de_genes), "unique genes to the STRING database..."))

# --- 3. QUERY STRING AND BUILD IGRAPH OBJECT ---
# Initialise the STRINGdb object for Gallus gallus (Taxonomy ID: 9031)
string_db <- STRINGdb$new(version = "12.0", 
                           species = 9031, 
                           score_threshold = 400, 
                           input_directory = "")

# Map our gene names to STRING's internal protein identifiers
mapped_genes <- string_db$map(de_genes, "gene_label", removeUnmappedRows = TRUE)

# Get the network as an igraph object, which is very stable
print("### Retrieving network as an igraph object...")
string_network_igraph <- string_db$get_graph()

# --- 4. ADD DATA TO THE NETWORK AND CUSTOMISE PLOT ---
# Get the logFC values for only the genes that are in our final network graph
# We match our data to the names of the vertices (nodes) in the graph
network_nodes <- V(string_network_igraph)$name
mapped_genes_in_network <- mapped_genes %>%
  dplyr::filter(STRING_id %in% network_nodes)

# Create a data frame that matches the graph's node order
node_data <- data.frame(STRING_id = V(string_network_igraph)$name) %>%
  dplyr::left_join(mapped_genes_in_network, by = "STRING_id")

# Set node colours based on logFC: red for up, blue for down
V(string_network_igraph)$color <- ifelse(node_data$logFC > 0, "red", "blue")

# --- 5. GENERATE AND SAVE THE NETWORK PLOT ---
print("### Generating and saving network plot with igraph...")
pdf("PPI_Network_Top_DE_Genes.pdf", width = 9, height = 9)
plot(string_network_igraph,
     vertex.size = 6,
     vertex.label = NA, # Hide labels on the plot to reduce clutter
     layout = layout_with_fr) # Use a nice force-directed layout
dev.off()


# --- 6. EXPORT FOR CYTOSCAPE (ADVANCED VISUALISATION) ---
print("### Exporting network data for Cytoscape...")
interaction_data <- string_db$get_interactions(mapped_genes$STRING_id)
readr::write_tsv(interaction_data, "cytoscape_edge_data.tsv")
readr::write_tsv(mapped_genes, "cytoscape_node_data.tsv")

print("### PPI network analysis complete. 'PPI_Network_Top_DE_Genes.pdf' has been saved.")
print("### For advanced visualisation, import 'cytoscape_edge_data.tsv' and 'cytoscape_node_data.tsv' into Cytoscape.")