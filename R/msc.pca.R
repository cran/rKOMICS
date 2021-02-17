#' Prinicple Component Analysis based on MSC
#'
#' The msc.pca allows you to perform Principle Component Analysis to summarize MSCs variation in all samples or in a subset of samples. 
#'
#' @usage msc.pca(clustmatrix, samples, groups, n = 20, labels = TRUE, title = NULL)
#' @param clustmatrix a cluster matrix
#' @param samples a vector containing the names of the samples. This can include all samples or it can be a subset. 
#' @param groups a vector specifying to which group (e.g. species) the samples belong to.
#' @param n number of clusters to select with highest contribution to PCA.
#' @param labels a logical parameter indicating whether to use labels on the PCA plot or not. Default is set to true.
#' @param title the title of the graph
#' @return
#'   \item{plot}{a PCA plot.}
#'   \item{eigenvalues}{a barplot showing percentage of explained variances.}
#'   \item{clustnames}{list of cluster names with highest contribution to PCA.}
#' @examples
#' data(matrices)
#' data(exData)
#' 
#' ### run function with all samples
#' res.pca <- lapply(matrices, function(x) msc.pca(x, samples = exData$samples, 
#'                   groups = exData$species, n=30, labels=FALSE, title=NULL))
#' 
#' res.pca$id93$eigenvalues
#' res.pca$id93$plot
#' 
#' ### use clusters with highest contribution to visualize in a heatmap
#' msc.heatmap(matrices[["id93"]][res.pca$id93$clustnames,], samples = exData$samples,
#'             groups = exData$species)
#'
#' ### run function with a subset of samples
#' ### you will be asked to confirm
#' table(exData$species)
#' hybrid <- which(exData$species=="hybrid")
#' # pca.subset <- msc.pca(clustmatrix = matrices[["id97"]], 
#' #                       samples = exData$samples[hybrid], 
#' #                       groups = exData$species[hybrid], labels = TRUE, 
#' #                       title = "PCA only with hybrids")
#' 
#' @importFrom factoextra fviz_pca_ind get_pca_var fviz_eig
#' @importFrom FactoMineR PCA
#' @import ggplot2
#' @importFrom utils menu
#' @export

msc.pca <- function(clustmatrix, samples, groups, n = 20, labels = TRUE, title = NULL) {

  
  
  ############# tests
  
  
  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (length(samples) != length(groups)) stop("ERROR: samples and groups are not of equal length")
  if (length(samples) != ncol(clustmatrix)) {
    answer <- utils::menu(c("Yes", "No"), 
                   title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
    if (answer!=1) stop("Function stopped")
  }
  
  
  ### extra toevoegen: PC kiezen

  pca <- FactoMineR::PCA(t(clustmatrix[,samples]),scale.unit=F, ncp=3, graph = F)

  
  
  if (labels == F) {
    plt <- factoextra::fviz_pca_ind(pca, geom.ind='point', col.ind = groups,
               addEllipses = F, legend.title="", alpha.ind = 0.8,
               pointsize = 4, invisible = "quali", title = title)
  } else {
    plt <- factoextra::fviz_pca_ind(pca, col.ind = groups, labelsize=3, 
                                    addEllipses = F, legend.title="", alpha.ind = 0.8,
                                    pointsize = 4, invisible = "quali",
                                    repel = TRUE, title = title)
  }

  contr <- factoextra::get_pca_var(pca)$contrib

  cl_index <- unique(order(contr[,1], decreasing = T)[1:n],
                         order(contr[,2], decreasing = T)[1:n],
                         order(contr[,3], decreasing = T)[1:n])
  cl_names <- rownames(clustmatrix[cl_index,samples])
  results <- list("plot" = plt,
                  "eigenvalues" = factoextra::fviz_eig(pca, addlabels = TRUE),
                  "clustnames" = cl_names)

  return(results)
}

