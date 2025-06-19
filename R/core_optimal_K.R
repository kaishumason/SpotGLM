
#' Compute Q-values from P-value List
#'
#' Applies Benjamini-Hochberg correction across a list of p-value vectors and reconstructs q-values in the same structure.
#'
#' @param pval_list A named list of numeric vectors containing p-values.
#'
#' @return A named list of numeric vectors with q-values corresponding to input p-values.
#'
#' @keywords internal
compute_qvals <- function(pval_list) {
  result <- vector("list", length(pval_list))
  names(result) <- names(pval_list)
  
  # Step 1: Flatten all p-values and apply BH correction
  all_pvals <- unlist(pval_list)
  all_qvals <- p.adjust(all_pvals, method = "BH")
  
  # Step 2: Reconstruct q-values list with original structure
  lengths_list <- sapply(pval_list, length)
  split_qvals <- split(all_qvals, rep(seq_along(lengths_list), lengths_list))
  
  # Step 3: Assign q-values to result list
  for (i in seq_along(split_qvals)) {
    result[[i]] <- split_qvals[[i]]
  }
  
  return(result)
}


#' Compute P-values for Pairwise Contrasts
#'
#' Computes Z-scores and p-values for contrasts between two covariates within a single cell type.
#'
#' @param input_list A list of model outputs with beta and covariance matrices.
#' @param cell_type Character. Cell type of interest.
#' @param effect_names Character vector of length 2 indicating the two covariates to contrast.
#' @param beta_name Name of the beta matrix in each model output (default: "beta_estimate").
#' @param covariance_name Name of the covariance matrix in each model output (default: "vcov").
#' @param sided Integer. 1 for one-sided, 2 for two-sided test.
#' @param direction Character. "pos" or "neg", used if `sided = 1`.
#'
#' @return A named list of p-values for each model in the input list.
#'
#' @keywords internal
compute_contrast <- function(input_list,cell_type,effect_names,beta_name = "beta_estimate", 
                             covariance_name = "vcov",
                             sided = 2, 
                             direction = "pos") {
  missing_counter = 0
  if(length(effect_names) != 2){
    stop("Effect names must be of length 2")
  }
  result = vector("list",length(input_list))
  test_statistics = vector("list",length(input_list))
  names(result) = names(input_list)
  names(test_statistics) = names(input_list)
  # Validate inputs
  if (!sided %in% c(1, 2)) stop("sided must be 1 or 2.")
  if (!direction %in% c("pos", "neg")) stop("direction must be 'pos' or 'neg'.")
  
  for (i in seq_along(input_list)) {
    beta_mat <- input_list[[i]][[beta_name]]
    vcov_mat <- input_list[[i]][[covariance_name]]
    
    
    #get cell type specific value
    cell_type_ind = which(colnames(beta_mat) == cell_type)
    effect_ind_1 = which(rownames(beta_mat) == effect_names[1])
    effect_ind_2 = which(rownames(beta_mat) == effect_names[2])
    #get number of rows and columns
    nX = nrow(beta_mat)
    nCT = ncol(beta_mat)
    #get vcov index for both effects
    effect_full_index_1 = (cell_type_ind - 1)*nX + effect_ind_1
    effect_full_index_2 = (cell_type_ind - 1)*nX + effect_ind_2
    
    
    if((length(cell_type_ind) > 1)){
      stop("Multiple cell type names match cell type")
    }else if( (length(effect_ind_1) > 1)){
      stop("Multiple feature names match effect name 1 ")
    }else if((length(effect_ind_2) > 1)){
      stop("Multiple feature names match effect name 2 ")
    }
    
    if((length(cell_type_ind) != 1) | (length(effect_ind_1) != 1) |(length(effect_ind_1) != 1) ){
      result[[i]] <- NA
      test_statistics[[i]] <- NA
      missing_counter = missing_counter + 1
      next
    }
    
    # Calculate Z-score
    difference = beta_mat[effect_ind_1,cell_type_ind] - beta_mat[effect_ind_2,cell_type_ind]
    standard_err = vcov_mat[effect_full_index_1,effect_full_index_1] 
    standard_err = standard_err + vcov_mat[effect_full_index_2,effect_full_index_2]
    standard_err = standard_err - 2*vcov_mat[effect_full_index_1,effect_full_index_2]
    standard_err = sqrt(standard_err)
    z <- difference / standard_err
    
    # Compute p-values
    if (sided == 2) {
      pvals <- 2 * pnorm(-abs(z))
    } else if (sided == 1) {
      if (direction == "pos") {
        pvals <- 1 - pnorm(z)
      } else {  # direction == "neg"
        pvals <- pnorm(z)
      }
    }
    # Save p-values to list
    result[[i]] <- pvals
    test_statistics[[i]] <- z
  }
  
  if(missing_counter>0){
    warning(paste0(missing_counter," entries were missing cell type or effect name"))
  }
  
  return(list(pvals = result,test_statistics = test_statistics))
}

#' Compute P-values for Single Covariate
#'
#' Computes Z-scores and p-values for a single covariate in a specific cell type.
#'
#' @param input_list A list of model outputs with beta and standard error matrices.
#' @param cell_type Character. Cell type of interest.
#' @param effect_name Character. Name of the covariate/effect of interest.
#' @param beta_name Name of the beta matrix in each model output (default: "beta_estimate").
#' @param standard_error_name Name of the standard error matrix (default: "standard_error_matrix").
#' @param sided Integer. 1 for one-sided, 2 for two-sided test.
#' @param direction Character. "pos" or "neg", used if `sided = 1`.
#'
#' @return A named list of p-values for each model in the input list.
#'
#' @keywords internal
compute_pvals <- function(input_list,cell_type,effect_name,beta_name = "beta_estimate", 
                          standard_error_name = "standard_error_matrix",
                          sided = 2, 
                          direction = "pos") {
  missing_counter = 0
  result = vector("list",length(input_list))
  test_statistics = vector("list",length(input_list))
  names(result) = names(input_list)
  names(test_statistics) = names(input_list)
  # Validate inputs
  if (!sided %in% c(1, 2)) stop("sided must be 1 or 2.")
  if (!direction %in% c("pos", "neg")) stop("direction must be 'pos' or 'neg'.")
  
  for (i in seq_along(input_list)) {
    beta_mat <- input_list[[i]][[beta_name]]
    se_mat <- input_list[[i]][[standard_error_name]]
    
    #get cell type specific value
    cell_type_ind = which(colnames(beta_mat) == cell_type)
    effect_ind = which(rownames(beta_mat) == effect_name)
    
    if((length(cell_type_ind) > 1)){
      stop("Multiple cell type names match cell type")
    }else if( (length(effect_ind) > 1) ){
      stop("Multiple feature names match effect name")
    }
    
    
    if((length(cell_type_ind)==0) | (length(effect_ind) == 0)){
      result[[i]] <- NA
      test_statistics[[i]] <- NA
      missing_counter = missing_counter + 1
      next
    }
    
    # Calculate Z-score
    z <- beta_mat[effect_ind,cell_type_ind] / se_mat[effect_ind,cell_type_ind]
    
    # Compute p-values
    if (sided == 2) {
      pvals <- 2 * pnorm(-abs(z))
    } else if (sided == 1) {
      if (direction == "pos") {
        pvals <- 1 - pnorm(z)
      } else {  # direction == "neg"
        pvals <- pnorm(z)
      }
    }
    # Save p-values to list
    result[[i]] <- pvals
    test_statistics[[i]] <- z
  }
  
  if(missing_counter>0){
    warning(paste0(missing_counter," entries were missing cell type or effect name"))
  }
  
  return(list(pvals = result,test_statistics = test_statistics))
}



#' Compute Significance for a Single Covariate
#'
#' Calculates p-values and adjusted q-values for a given effect in a specific cell type.
#'
#' @param input_list A list of model outputs with beta and standard error matrices.
#' @param cell_type Character. Cell type of interest.
#' @param effect_name Character. Name of the covariate/effect of interest.
#' @param beta_name Name of the beta matrix (default: "beta_estimate").
#' @param standard_error_name Name of the standard error matrix (default: "standard_error_matrix").
#' @param sided Integer. 1 for one-sided, 2 for two-sided test.
#' @param direction Character. "pos" or "neg", used if `sided = 1`.
#'
#' @return A data frame with columns: name, pval, and qval.
#'
#' @export
compute_significance <- function(input_list,cell_type,effect_name,beta_name = "beta_estimate", 
                                 standard_error_name = "standard_error_matrix", 
                                 sided = 2, 
                                 direction = "pos") {
  if(is.null(names(input_list))){
    names(input_list) = paste0("Test_",c(1:length(input_list)))
  }
  
  #Step 1: Compute pvalues
  pvals = compute_pvals(input_list = input_list,cell_type = cell_type,effect_name = effect_name, beta_name = beta_name,
                        standard_error_name = standard_error_name,sided = sided,
                        direction = direction)
  #Step 2: Compute qvalues
  qvals = compute_qvals(pval_list = pvals$pvals)
  
  #return matrix of results 
  result_mat <- do.call(rbind, lapply(seq_along(pvals$pvals), function(i) {
    data.frame(
      name = names(pvals$pvals)[i],
      cell_type = cell_type,
      effect = effect_name,
      sided = sided,
      direction = direction,
      test_statistic = pvals$test_statistics[[i]],
      pval = pvals$pvals[[i]],
      qval = qvals[[i]],
      stringsAsFactors = FALSE
    )
  }))
  
  # Return the matrix (as a data frame)
  return(result_mat)
  
}

#' Compute Significance for Pairwise Covariate Contrast
#'
#' Calculates p-values and adjusted q-values for contrasts between two effects within a specific cell type.
#'
#' @param input_list A list of model outputs with beta and covariance matrices.
#' @param cell_type Character. Cell type of interest.
#' @param effect_names Character vector of length 2 indicating the covariates to contrast.
#' @param beta_name Name of the beta matrix (default: "beta_estimate").
#' @param covariance_name Name of the covariance matrix (default: "vcov").
#' @param sided Integer. 1 for one-sided, 2 for two-sided test.
#' @param direction Character. "pos" or "neg", used if `sided = 1`.
#'
#' @return A data frame with columns: name, pval, and qval.
#'
#' @export
compute_contrast_significance = function(input_list,cell_type,effect_names,beta_name = "beta_estimate", 
                                         covariance_name = "vcov", 
                                         sided = 2, 
                                         direction = "pos") {
  if(is.null(names(input_list))){
    names(input_list) = paste0("Test_",c(1:length(input_list)))
  }
  
  #Step 1: Compute pvalues
  pvals = compute_contrast(input_list = input_list,cell_type = cell_type,effect_names = effect_names,beta_name = beta_name,
                           covariance_name = covariance_name,sided = sided,
                        direction = direction)
  #Step 2: Compute qvalues
  qvals = compute_qvals(pval_list = pvals$pvals)
  
  #return matrix of results 
  result_mat <- do.call(rbind, lapply(seq_along(pvals$pvals), function(i) {
    data.frame(
      name = names(pvals$pvals)[i],
      cell_type = cell_type,
      effect_1 = effect_names[1],
      effect_2 = effect_names[2],
      sided = sided,
      direction = direction,
      test_statistic = pvals$test_statistics[[i]],
      pval = pvals$pvals[[i]],
      qval = qvals[[i]],
      stringsAsFactors = FALSE
    )
  }))
  
  # Return the matrix (as a data frame)
  return(result_mat)
  
}




#' Read MERFISH Example Data from GitHub
#'
#' Downloads and loads MERFISH example data from a public GitHub repository.
#' The function retrieves cell type annotations, spatial region labels, and
#' gene expression count matrices (in up to 10 chunks), then returns them as a list.
#'
#' @param num_chunks Integer between 0 and 10. Controls how many count matrix chunks
#' to download and load. Defaults to 10 (all chunks). Use 0 to skip loading counts.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{counts}{A matrix or data frame of gene expression counts (combined from multiple chunks).}
#'   \item{regions}{A vector or factor of region labels for each cell.}
#'   \item{CT}{A vector or factor of cell types for each cell.}
#' }
#'
#' @details This function requires an internet connection to download data from
#' \url{https://github.com/kaishumason/SpotGLM-Example-Data}. The function assumes the
#' repository structure and filenames are consistent with the expected format.
#'
#' @examples
#' \dontrun{
#' data_list <- read_merfish_data(num_chunks = 5)
#' head(data_list$counts)
#' table(data_list$CT)
#' }
#'
#' @export
read_merfish = function(num_chunks = 10){
  #get cell types
  # Define raw GitHub URL
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/merfish/cell_types.rds"
  
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  cell_types <- readRDS(temp_file)
  
  
  #get regions
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/merfish/regions.rds"
  
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  regions <- readRDS(temp_file)
  
  #get counts
  n = min(as.integer(num_chunks),10)
  if(n<0){
    stop("number of chunks should be an integer between 0(no counts) and 10 (all counts)")
  }
  data = vector("list",n)
  for(j in c(1:n)){
    url <- paste0("https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/merfish/counts_",j,".rds")
    # Temporary file to store the .rds
    temp_file <- tempfile(fileext = ".rds")
    
    # Download the file (use mode = "wb" for binary)
    download.file(url, destfile = temp_file, mode = "wb")
    
    data[[j]] = readRDS(temp_file)
  }
  data = do.call(rbind, data)
  
  return(list(counts = data,regions = regions, CT = cell_types))
}



#' Read Visium HD Example Data from GitHub
#'
#' Downloads and loads spatial transcriptomics data from the Visium HD example
#' dataset hosted on a public GitHub repository. This includes spatial coordinates,
#' gene expression counts, deconvolution results, and effective niche estimates.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{coords}{A matrix or data frame of spatial coordinates (x and y) for each spot.}
#'   \item{niche}{A matrix, data frame, or list representing effective niche composition per spot.}
#'   \item{deconv}{A matrix or data frame of cell type deconvolution proportions per spot.}
#'   \item{counts}{A gene expression count matrix (genes × spots).}
#' }
#'
#' @details This function requires an internet connection to download data from
#' \url{https://github.com/kaishumason/SpotGLM-Example-Data}. The repository must contain
#' the files \code{coords.rds}, \code{deconv_matrix.rds}, \code{count_matrix_subset.rds},
#' and \code{niche.rds} in the \code{visiumHD} folder.
#'
#' @examples
#' \dontrun{
#' data_list <- read_visiumHD()
#' head(data_list$coords)
#' dim(data_list$counts)
#' }
#'
#' @export
read_visiumHD = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/coords.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file)[,4:5])
  
  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/deconv.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file))
  
  
  #read counts
  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/counts.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  counts <- readRDS(temp_file)
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/niche.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  niche <- readRDS(temp_file)
  
  #read library size
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/library_size.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  library_size <- readRDS(temp_file)
  
  
  return(list(coords = coords,niche = niche, deconv = deconv, counts = counts,library_size = library_size))
  
}






#' Read Visium Example Data from GitHub
#'
#' Downloads and loads Visium spatial transcriptomics data from a GitHub repository.
#' Includes coordinates, deconvolution, effective niche covariates, library sizes, and gene counts.
#'
#' @return A list containing:
#' \describe{
#'   \item{coords}{Matrix of spatial coordinates.}
#'   \item{niche}{Effective niche covariate matrix.}
#'   \item{deconv}{Cell type deconvolution matrix.}
#'   \item{counts}{Gene expression count matrix.}
#'   \item{library_size}{Vector of library sizes per spot.}
#' }
#'
#' @examples
#' \dontrun{
#' visium_data <- read_example_visium_colorectal_cancer_data()
#' str(visium_data)
#' }
#'
#' @export
read_example_visium_colorectal_cancer_data = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/coords.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file))
  
  
  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/deconv.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file))
  
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/effective_niche_covariates.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  niche <- readRDS(temp_file)
  
  
  
  #read library size
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/library_size.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  library_size <- readRDS(temp_file)
  
  
  #get counts
  data = vector("list",4)
  for(j in c(1:4)){
    url <- paste0("https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/counts_",j,".rds")
    # Temporary file to store the .rds
    temp_file <- tempfile(fileext = ".rds")
    
    # Download the file (use mode = "wb" for binary)
    download.file(url, destfile = temp_file, mode = "wb")
    
    data[[j]] = readRDS(temp_file)
  }
  data = do.call(cbind, data)
  
  
  
  return(list(coords = coords,niche = niche, deconv = deconv, counts = data,library_size = library_size))
  
}



#' Read Spatial ATAC-Seq Example Data
#'
#' Loads spatial ATAC-seq data from a GitHub repository, including coordinates, deconvolution,
#' region-level features, and motif score matrices.
#'
#' @return A list containing:
#' \describe{
#'   \item{coords}{Spatial coordinates matrix.}
#'   \item{regions}{Matrix or data frame of spatial regions per spot.}
#'   \item{deconv}{Cell type deconvolution matrix.}
#'   \item{motif_scores}{Matrix of motif activity scores per region.}
#' }
#'
#' @examples
#' \dontrun{
#' atac_data <- read_example_spatial_mouse_brain_atac_data()
#' names(atac_data)
#' }
#'
#' @export
read_example_spatial_mouse_brain_atac_data = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/coord.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file))
  
  
  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/deconvolution.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file)$mat)
  
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/region_matrix.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  regions <- readRDS(temp_file)
  
  
  
  
  #get counts
  data = vector("list",3)
  for(j in c(1:3)){
    url <- paste0("https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/scores_",j,".rds")
    # Temporary file to store the .rds
    temp_file <- tempfile(fileext = ".rds")
    
    # Download the file (use mode = "wb" for binary)
    download.file(url, destfile = temp_file, mode = "wb")
    
    data[[j]] = readRDS(temp_file)
  }
  data = do.call(cbind, data)
  
  
  
  return(list(coords = coords,regions = regions, deconv = deconv, motif_scores = data))
  
}



#' Read Spatial Long-Read RNA-Seq Data
#'
#' Loads spatial long-read RNA-seq data from our public GitHub repository. Includes coordinates,
#' region annotations, deconvolution, library sizes, and expression matrices for genes and isoforms.
#'
#' @return A list containing:
#' \describe{
#'   \item{coords}{Matrix of spatial coordinates.}
#'   \item{regions}{Spatial region annotations.}
#'   \item{deconv}{Cell type deconvolution matrix.}
#'   \item{library_size}{Vector of library sizes per spot.}
#'   \item{total_gene_expression}{Matrix of total gene expression.}
#'   \item{isoform_expression}{Matrix of isoform-level expression.}
#' }
#'
#' @examples
#' \dontrun{
#' long_read_data <- read_example_spatial_olfactory_bulb_long_read_data()
#' str(long_read_data)
#' }
#'
#' @export
read_example_spatial_olfactory_bulb_long_read_data = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/coords.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file))
  
  
  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/deconvolution.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file))
  
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/region_matrix.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  regions <- readRDS(temp_file)
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/library_size.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  library_size <- readRDS(temp_file)
  
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/total_gene_expression.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  total_gene_expression <- readRDS(temp_file)
  
  
  
  
  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/isoform_expression.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")
  
  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")
  
  # Read the RDS file
  isoform_expression <- readRDS(temp_file)
  
  
  
  
  
  
  return(list(coords = coords,regions = regions, deconv = deconv,library_size = library_size,
              total_gene_expression = total_gene_expression ,isoform_expression = isoform_expression))
  
}





