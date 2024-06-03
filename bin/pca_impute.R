#!/usr/bin/env Rscript

## Imputation of missing data using PCA
box::use(missMDA[...])



## Main funciton
#' @export
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
            datt[, !(colnames(datt) %in% columns), drop = FALSE], 
            imputed$completeObs)
      } else {
      	## Return completed dataset
        res <- imputed$completeObs
      }
    }
  }
    return(as.data.frame(res))
}




