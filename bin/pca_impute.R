#!/usr/bin/env Rscript

## Imputation of missing data using PCA

library(missMDA)
library(optparse)
library(readxl)


selected_meta_columns = c("MeanTemperature_Avg31day",
                    "PrecipitationSum_Avg31day",
                    "pH_H2O",
                    "Electrical_conductivity",
                    "Carbonate_content",
                    "Phosphorus_content",
                    "CN",
                    "Clay_content_imputed",
                    "Organic_carbon_imputed",
                    "H2O_content_volumetric_imputed",
                    "Annual_Precipitation",
                    "Annual_Mean_Temperature",
                    "Bulk_Density_0_10_cm_imputed",
                    "Bulk_Density_10_20_cm_imputed")

selected_string_columns = c("SampleID","LC_simplest")
                    

# Define the option parser
option_list <- list(
    make_option(c("--metadata_file"),
        type = "character", default = "/home/minigonche/Dropbox/Projects/TartuU/community_network_analyzer/LUCAS_Funct/metadata_for_samples_2.xlsx", # nolint
        help = "Location of the metadata file"
    )
)

# Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))


metadata_file <- opt$metadata_file

# Reads
df_meta <- read_excel(metadata_file)

df_meta = as.data.frame(df_meta)


## Main function

## Main funciton
pca_impute <- function(datt, columns = NULL){
  # columns = vector of column names used to construct PCA

  ## Subset data for PCA
  if(is.null(columns)){
    x <- datt
  } else {
    x <- datt[, columns]
  }

  ## Estimate number of missing values per row
  nmissing <- function(x) sum(is.na(x))
  row.na <- apply(x, 1, nmissing)

  ## If there's no missing values, just return the data
  if(sum(row.na) == 0){ return(datt) }

  ## Otherwise, perform PCA-based missing data imputation
  if(sum(row.na) > 0){
  
  ## The number of components which leads to the smallest MSEP
  nb <- try(estim_ncpPCA(x,
  	      scale     = TRUE,
          ncp.min   = 1,
          ncp.max   = floor(ncol(x)*0.8),
          method    = "Regularized", # regularized iterative PCA algorithm for data imputation
          method.cv = "gcv",         # generalised cross-validation
          nbsim = 1000), silent=TRUE)
          
    if(class(nb) == "try-error") {
      cat("....Imputatuion failed\n")
      res <- datt
    }
    if(class(nb) != "try-error"){
      imputed <- imputePCA(x,
      	     scale = TRUE,
              ncp = nb$ncp,
              method = "Regularized",
              nb.init = 10,           # number of random initializations        
              row.w = 1/(row.na+1))   # row weights proportional to the number of missing vals

      if(!is.null(columns)){
        ## Recover columns that were not used in the analysis
        res <- cbind(
            datt[, !(colnames(datt) %in% columns)], 
            imputed$completeObs)
      } else {
      	## Return completed dataset
        res <- imputed$completeObs
      }
    }
  }
    return(as.data.frame(res))
}


# Function Call
df_meta_new = pca_impute(df_meta, selected_meta_columns)

# reattaches the necessary columns
df_meta_new = cbind(df_meta[,selected_string_columns], df_meta_new)

# Export
write.csv(df_meta_new, gsub(".xlsx", "_imputed.csv", input_file), row.names = FALSE)

