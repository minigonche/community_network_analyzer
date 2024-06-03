#!/usr/bin/env Rscript

## Usage
# Rscript network_builder.R \
#  --input_file gene_abundances_long.csv \
#  --edge_sign both \
#  --metadata_file sample_metadata.xlsx \
#  --replace_missing FALSE \
#  --metadata_cols "CN,H2O_content_volumetric,Annual_Precipitation,Annual_Mean_Temperature"


load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("SPRING")
load_pckg("tidyr")
load_pckg("dplyr")
load_pckg("vegan")
load_pckg("igraph")
load_pckg("optparse")
load_pckg("readxl")
load_pckg("missMDA")


box::use(./pca_impute[...])


seed <- 42

# Define the option parser
option_list <- list(
    make_option(c("--input_file"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/work/ba/10afc8be37ef894256ad78789ff100/P_cycle-grassland.csv",
        help = "Location of the input file"
    ),
    make_option(c("--only_positive"),
        type = "logical", default = FALSE, 
        help = "If the graph should only include positive weights"
    ),  
    make_option(
        c("--edge_sign"), 
        type = "character", 
        default = "both", 
        help = "Mode of operation: 'positive', 'negative', or 'both'", 
        #choices = c("positive", "negative", "both"),
    ),
    make_option(c("--metadata_file"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/LUCAS_Funct/metadata_for_samples_2.xlsx", # nolint
        help = "Location of the metadata file"
    ),
    make_option(c("--metadata_cols"),
        type = "character",
        default = "MeanTemperature_Avg31day,PrecipitationSum_Avg31day,pH_H2O,Electrical_conductivity,Carbonate_content,Phosphorus_content,CN,Clay_content_imputed,Organic_carbon_imputed,H2O_content_volumetric_imputed,Annual_Precipitation,Annual_Mean_Temperature,Bulk_Density_0_10_cm_imputed,Bulk_Density_10_20_cm_imputed",
        help = "Selected metadata columns (comma-separated)"
    ),
    make_option(c("--replace_missing"),
        type = "logical", default = TRUE, 
        help = "Replace missing metadata"
    )
)

# Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Starts the parameters
input_file <- opt$input_file
edge_sign <- opt$edge_sign
metadata_file <- opt$metadata_file
selected_meta_columns <- opt$metadata_cols
replace_missing <- opt$replace_missing

# Load metadata
if(grepl("\\.xlsx$", metadata_file) || grepl("\\.xls$", metadata_file)){
    df_meta <- read_excel(metadata_file)
} else {
    df_meta <- read.csv(metadata_file)
}



# Split metadata columns
if(selected_meta_columns %in% "NA"){ selected_meta_columns <- NA }
if(is.na(selected_meta_columns)){
    selected_meta_columns <- colnames(df_meta)[ ! colnames(df_meta) %in% "SampleID" ]
} else {
    selected_meta_columns <- unique(strsplit(selected_meta_columns, split = ",")[[1]])
}

# Filters out columns not in meta file
missing_cols = setdiff(selected_meta_columns, colnames(df_meta))
selected_meta_columns = intersect(selected_meta_columns, colnames(df_meta))
if(length(missing_cols) > 0){
    cat("\nWARNING: Some metadata columns are not included in the metadata file provided:\n")
    cat("  ", paste(missing_cols, collapse = ", "), "\n")
    cat(".. They will not be included in analysis\n")
}

# Subset metadata
df_meta <- df_meta  %>%
    select("SampleID", all_of(selected_meta_columns))

df_meta <- as.data.frame(df_meta)


## Handle missing data in metadata
nacolz <- colSums(is.na(df_meta)) != 0
if(any(nacolz)){
    cat("\nWARNING: There are ", sum(nacolz), "columns with missing data:\n")
    cat("  ", paste(colnames(df_meta)[ nacolz ], collapse = ", "), "\n")


    if(replace_missing == FALSE){
            ## Removes NA Columns
            cat(".. These columns will be removed\n")
            df_meta <- df_meta[, ! nacolz ]
    } else {
        # Replace missing values in metadata
            cat(".. Imputing missing values\n")
            # Function Call
            df_meta <- pca_impute(df_meta, selected_meta_columns)

    }
}

# Metadata variables that will be used for ordination
final_meta_columns <- colnames(df_meta)
final_meta_columns <- final_meta_columns[ ! final_meta_columns %in% "SampleID" ]

# Load gene abundances, Reads and removes unused columns
df <- read.csv(input_file) %>%
    select("SampleID", "Gene", "Abundance")


# Validate sample sets
nonoverlaping_samples <- setdiff(unique(df$SampleID), unique(df_meta$SampleID))
if(length(nonoverlaping_samples) > 0){
    cat("\nWARNING: some samples are missing from metadata or gene-abundance table:\n")
    cat("  n = ", length(nonoverlaping_samples), "; ", paste(nonoverlaping_samples, collapse = ", "), "\n")
    cat(".. Subsetting to a common sample set\n")

    samples_in_common <- intersect(unique(df$SampleID), unique(df_meta$SampleID))
    if(length(samples_in_common) == 0){ stop("ERROR: no samples in common!\n") }
    df <- df[ df$SampleID %in% samples_in_common,  ]
    df_meta <- df_meta[ df_meta$SampleID %in% samples_in_common,  ]
}

rownames(df_meta) <- df_meta$SampleID

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
col_sums <- col_sums[col_sums == 0]
all_zero_genes <- names(col_sums)

df <- df %>%
  select(-all_of(all_zero_genes))

# Converts to matrix 
X <- as.matrix(df)

# Runs Spring
cat("\nRunning SPRING\n")
fit.spring <- suppressPackageStartupMessages(
    SPRING(
    X,
    Rmethod = "approx",
    quantitative = TRUE,
    lambdaseq = "data-specific",
    nlambda = 50,
    rep.num = 50,
    verbose = FALSE,
    seed = seed
    ))

# Uses SpiecEasi to invert the similarity matrix
cat("Inverting similarity matrix\n")
opt.K <- fit.spring$output$stars$opt.index
pcor.K <- as.matrix(
    SpiecEasi::symBeta(
    fit.spring$output$est$beta[[opt.K]],
    mode = "maxabs"
    )
)

# Creates Object
df_matrix <- pcor.K

if(edge_sign == "positive")
    df_matrix[df_matrix < 0] <- 0
    
if(edge_sign == "negative")
    df_matrix[df_matrix > 0] <- 0

## Assign column and row names
colnames(df_matrix) <- colnames(df)
rownames(df_matrix) <- colnames(df)

# Creates the Graph
cat("Costructing a graph\n")
g <- graph_from_adjacency_matrix(df_matrix, weighted = TRUE, mode = "undirected")

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
cat("Estimating modularity\n")
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

cat("Running db-RDA\n")
cap <- capscale(
  formula,
  distance = "euclidean",
  data = df_meta)

cat(".. number of constrained axes: ", length(cap$CCA$eig), "\n")

## Extract scores
capscores <- data.frame(
  GeneID = rownames(cap$CCA$v),
  scores(cap, display = "species", choices = 1:length(cap$CCA$eig)))

## Add CAP scores as node attributes
cat("Adding db-RDA scores to graph\n")
node_labs <- V(g)$name
for(i in seq(1, length(cap$CCA$eig))){
    capCol = paste("CAP", i, sep = "")
    
    # Adds value
    g <- set_vertex_attr(
        g,
        capCol,
        index = node_labs, capscores[ match(x = node_labs, table = capscores$GeneID), capCol ]
        )

    # Adds Assortativity
    g <- set_graph_attr(
        g,
        paste("assortativity", capCol, sep = "_"),
        igraph::assortativity(pos_g, values = vertex_attr(g, capCol, index = V(g)))
        )

    # Add the biplot values for projection of the CAP values
    for(col in rownames(cap$CCA$biplot)){
        g <- set_graph_attr(g, paste(capCol, col, sep = "-"), cap$CCA$biplot[col, capCol])
    }

}

# Writes Graph
cat("\nExporting graph\n")
write_graph(g,  gsub(".csv", ".graphml", input_file), "graphml")

# Write auxiliary data in R format
save(g, cap, fit.spring,
    file = gsub(".csv", ".RData", input_file),
    compress = "xz")
