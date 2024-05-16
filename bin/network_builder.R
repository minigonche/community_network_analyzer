#!/usr/bin/env Rscript

# SPRING Testing
library(SPRING)
library(tidyr)
library(dplyr)
library(vegan)
library(igraph)
library(optparse)
library(readxl)
library(missMDA)


seed <- 42

# Define the option parser
option_list <- list(
    make_option(c("--input_file"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/work/3f/cb397522d28413b7c04ca0703fe757/S_cycle-grassland.csv",
        help = "Location of the input file"
    ),
    make_option(c("--only_positive"),
        type = "logical", default = FALSE, 
        help = "If the graph should only include positive weights"
    ),  
    make_option(c("--metadata_file"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/LUCAS_Funct/metadata_for_samples.xlsx", # nolint
        help = "Location of the metadata file"
    )  
)

# Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Starts the parameters
input_file <- opt$input_file
only_positive <- opt$only_positive
metadata_file <- opt$metadata_file


# Reads Metadata
df_meta <- read_excel(metadata_file)  %>%
    select("SampleID", all_of(selected_meta_columns))

df_meta <- as.data.frame(df_meta)

rownames(df_meta) = df_meta$SampleID

# Removes NA Columns
df_meta <- df_meta[, colSums(is.na(df_meta)) == 0]

final_meta_columns = colnames(df_meta)
final_meta_columns = final_meta_columns[final_meta_columns != "SampleID"]


# Reads and removes unused columns
df <- read.csv(input_file) %>%
    select("SampleID", "Gene","Abundance")

# Pivots
# Genes as Columns
df <- df %>%
    pivot_wider(names_from = Gene, values_from = Abundance) 
    
# Aligns Meta and Abundance
df_meta <- df_meta[df_meta$SampleID %in% df$SampleID, ]
df_meta <- df_meta[match(df_meta$SampleID, df$SampleID), ]

# Drops SampleID from abundance 
df <- df %>% select(-SampleID)

# Drops genes that are not present in any sample
col_sums <- colSums(df)
col_sums = col_sums[col_sums == 0]
all_zero_genes = names(col_sums)

df <- df %>%
  select(-all_of(all_zero_genes))

# Converts to matrix 
X <- df %>% as.matrix()

# Runs Spring
fit.spring <- SPRING(
    X,
    Rmethod = "approx",
    quantitative = TRUE,
    lambdaseq = "data-specific",
    nlambda = 50,
    rep.num = 50,
    verbose = FALSE,
    seed = seed
)

# Uses SpiecEasi to invert the similarity matrix
opt.K <- fit.spring$output$stars$opt.index
pcor.K <- as.matrix(
    SpiecEasi::symBeta(
    fit.spring$output$est$beta[[opt.K]],
    mode = "maxabs"
    )
)

# Creates Object
df_matrix <- pcor.K

if (only_positive)
    df_matrix[df_matrix < 0] <- 0

# Converts to Tibble 
df_matrix <- as_tibble(df_matrix)

colnames(df_matrix) <- colnames(df)
rownames(df_matrix) <- colnames(df)

# Creates the Graph
g <- graph_from_adjacency_matrix(as.matrix(df_matrix), weighted = TRUE, mode = "undirected")

# Creates the different versions of the graph
# Only Positive
# All statistics are extracted using only positive edges
pos_g <- delete_edges(g, E(g)[E(g)$weight < 0])

# Extracts Statistics
# Abundance
V(g)$max_abundance <- apply(df, 2, max)
V(g)$min_abundance <- apply(df, 2, min)

# General Graph Statistics
g = set_graph_attr(g, "density", edge_density(pos_g))
g = set_graph_attr(g, "avg_weight", mean(E(pos_g)$weight))

# Modularity
communities <- igraph::cluster_fast_greedy(pos_g, 
                                           weights = E(pos_g)$weight, 
                                           merges = TRUE, 
                                           modularity = TRUE, 
                                           membership = TRUE)

g = set_graph_attr(g, "modularity", modularity(communities))
g = set_graph_attr(g, "num_clusters", max(membership(communities)))
V(g)$membership = membership(communities)

# Assortativity

## mCLR
mdf <- SPRING::mclr(df, base = exp(1), tol = 1e-16, eps = NULL, atleast = 1)

## Constrained ordination (db-RDA)
formula <- as.formula(paste("mdf", "~", paste(final_meta_columns, collapse = " + ")))

cap <- capscale(
  formula,
  distance = "euclidean",
  data = df_meta)

## Extract scores
capscores <- data.frame(
  GeneID = rownames(cap$CCA$v),
  scores(cap, display = "species", choices = seq(1, length(final_meta_columns))))

## Add CAP scores as node attributes
node_labs <- V(g)$name
for(i in seq(1, length(final_meta_columns)))
{
    capCol = paste("CAP", i, sep = "")
    # Adds value
    g <- set_vertex_attr(g, capCol, index = node_labs, capscores[ match(x = node_labs, table = capscores$GeneID), capCol ])

    # Adds Assortativity
    g <- set_graph_attr(g, paste("assortativity",capCol, sep = "_"), igraph::assortativity(pos_g, values = vertex_attr(g, capCol, index = V(g))))

    # Add the biplot values for projection of the CAP values
    for(col in final_meta_columns)
    {
        g <- set_graph_attr(g, paste(capCol, col, sep = "-"), cap$CCA$biplot[col, capCol])
   
    }

}


# Writes Graph
write_graph(g,  gsub(".csv", ".graphml", input_file), "graphml")

