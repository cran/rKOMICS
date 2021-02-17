#' Specific MSC
#'
#' The msc.subset allows you to find specific MSC for a certain subset of samples.
#'
#' @usage msc.subset(clustmatrix, subset)
#' @param clustmatrix a cluster matrix
#' @param subset a numerical vector indicating which subset of samples to include.
#' @return
#'   \item{clustnumbers}{a vector containing the specific MSC names.}
#'   \item{freq}{frequency values of those specific MSC in the subset of samples.}
#'   \item{matrix}{a subset of the cluster matrix containing only those specific MSC. All samples, not in the subset, should have a value of 0 meaning the MSC is absent.}
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


  ############# tests


  if (!is.matrix(clustmatrix)) stop("ERROR: clustmatrix should be a matrix")
  if (!is.numeric(subset)) stop("ERROR: subset should be a numerical vector")


  ############# spec MSC


  spec_MSC <- list()
  spec_MSC$sum <- sum(apply(clustmatrix[,-subset], 1, function(x) sum(x)) == 0)
  spec_MSC$clustnumbers <- names(which(apply(clustmatrix[,-subset], 1, function(x) sum(x)) == 0))
  index <- which(apply(clustmatrix[,-subset], 1, function(x) sum(x)) == 0)
  spec_MSC$matrix <- clustmatrix[index,]
  spec_MSC$freq <- rowSums(clustmatrix[index,,drop=F])

  return(spec_MSC)

}

