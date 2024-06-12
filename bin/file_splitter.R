#!/usr/bin/env Rscript

# Script for splitting the given file into the different land covers
library(readxl)
library(tidyr)
library(dplyr)
library(optparse)



# Command
# Rscript file_splitter.R

# Define the option parser
option_list <- list(
    make_option(c("--cycle"),
        type = "character", default = "S_cycle",
        help = "Name of the cycle"
    ),
    make_option(c("--input_folder"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/LUCAS_Funct", # nolint
        help = "Location of input folder"
    ),
    make_option(c("--level"),
        type = "character", default = "L1",
        help = "Level to analyze: L1, L2, L3 or L4"
    ),
    make_option(c("--metadata_file"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/LUCAS_Funct/metadata_for_samples.xlsx", # nolint
        help = "Location of the metadata file"
    ),
    make_option(c("--split_by_lc"),
        type = "logical", default = TRUE, 
        help = "Split the file by land cover"
    )
)

# Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))


input_cycle <- opt$cycle
input_folder <- opt$input_folder
level <- opt$level
metadata_file <- opt$metadata_file
split_by_lc <- opt$split_by_lc


# Start Scripts
input_file <- file.path(input_folder, paste(input_cycle, ".xlsx", sep = ""))

land_cover_column <- "LC_simplest"
land_cover_options <- c(
    "cropland",
    "grassland",
    "woodland_coniferous",
    "woodland_deciduous"
)

# Reads Metadata
if(grepl("\\.xlsx$", metadata_file) || grepl("\\.xls$", metadata_file)){
    df_meta <- read_excel(metadata_file)
} else {
   df_meta <- read.csv(metadata_file)
}


df_meta <- df_meta %>% select("SampleID", all_of(land_cover_column))

# Reads the file
df <- read_excel(input_file)
# Filters out non relevant info
df <- df %>%
    filter(AggregLevel == level) %>%
    select("Gene", "SampleID", "Abundance")
# Merges
df <- df %>% left_join(df_meta, by <- join_by(SampleID == SampleID))

# Outputs all all land covers
df_temp <- df %>%
        select("Gene", "SampleID", "Abundance", all_of(land_cover_column))
# Write the file
write.csv(df_temp, paste(input_cycle, "-", "all_lc", ".csv", sep = ""), row.names = FALSE)

if(split_by_lc)
{
    # Filters by land cover
    for (lc in land_cover_options) {
        df_temp <- df %>%
            filter(get({{ land_cover_column }}) == lc) %>%
            select("Gene", "SampleID", "Abundance", all_of(land_cover_column))
        # Write the file
        write.csv(df_temp, paste(input_cycle, "-", lc, ".csv", sep = ""), row.names = FALSE)
    }

}

