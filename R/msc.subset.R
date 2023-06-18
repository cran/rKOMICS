#' Specific Minicircle Sequence Classes
#'
#' The msc.subset function allows you to identify specific Minicircle Sequence Classes (MSCs) for a subset of samples based on the output of the msc.matrix function. It helps in extracting and analyzing MSCs that are present in a particular subset of samples.
#'
#' @usage msc.subset(clustmatrix, subset)
#' @param clustmatrix a cluster matrix obtained from the msc.matrix function. The cluster matrix represents the presence or absence of MSCs in each sample.
#' @param subset a numerical vector indicating the subset of samples for which you want to identify specific MSCs. The values in the subset vector correspond to the indices of the samples to be included.
#' @return
#'   \item{clustnumbers}{a vector containing the names of the specific MSCs present in the subset of samples. These are the MSCs that are found in the indicated subset.}
#'   \item{freq}{frequency values indicating the occurrence of the specific MSCs in the subset of samples. These values represent the number of times each MSC appears in the subset.}
#'   \item{matrix}{a subset of the cluster matrix containing only the specific MSCs found in the subset of samples. For samples not included in the subset, the values in the matrix should have the value 0, indicating the absence of the MSC.}
#'   \item{sum}{the total number of MSC found in the indicated subset of samples.}
#' @examples
#' data(matrices)
#' data(exData)
#'
#' ### selecting a group of samples e.g. all L. peruviana species
#' Lpe <- which(exData$species == "L. peruviana")
#' 
#' ### run function
#' specific <- msc.subset(matrices[["id97"]], subset = Lpe)
#' 
#' ### visualize results (check if it is indeed specific)
#' heatmap(specific$matrix) # or:
#' msc.heatmap(specific$matrix, samples = exData$samples, groups = exData$species)
#'
#' ### find specific MSC with highest frequency
#' which.max(specific$freq)
#'
#' @export

msc.subset <- function(clustmatrix, subset) {


  #############   0   Tests   #############


  if (!is.matrix(clustmatrix)) stop("ERROR: clustmatrix should be a matrix")
  if (!is.numeric(subset)) stop("ERROR: subset should be a numerical vector")


  #############   1   Specific MSC   #############


  spec_MSC <- list()
  index <- apply(clustmatrix[, -subset], 1, function(x) sum(x) == 0)
  spec_MSC$sum <- sum(index)
  spec_MSC$clustnumbers <- names(which(index))
  spec_MSC$matrix <- clustmatrix[index, ]
  spec_MSC$freq <- rowSums(clustmatrix[index, , drop = FALSE])

  
  #############   2   Return specific MSC   #############
  
  
  return(spec_MSC)

}

