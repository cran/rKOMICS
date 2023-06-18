#' Prinicple Component Analysis based on MSC
#'
#' The msc.pca function allows you to perform Principle Component Analysis (PCA) to summarize the variation of Minicircle Sequence Classes (MSCs) in all samples or in a subset of samples.
#'
#' @usage msc.pca(clustmatrix, samples, groups, n = 20, labels = TRUE, title = NULL)
#' @param clustmatrix a cluster matrix obtained from the msc.matrix function. The cluster matrix represents the presence or absence of MSCs in each sample, where rows represent MSCs and columns represent samples.
#' @param samples a vector containing the names of the samples. This can include all samples or a subset of samples that you want to analyze.
#' @param groups a vector specifying the groups (e.g., species) to which the samples belong.
#' @param n the number of clusters to select with the highest contribution to PCA. By default, it is set to 20.
#' @param labels a  logical parameter indicating whether to use labels on the PCA plot or not. If set to TRUE (default), the plot will display sample labels.
#' @param title the  title of the graph. You can provide a title for the PCA plot if desired.
#' @return
#'   \item{plot}{a PCA plot that visualizes the clustering of samples based on the presence/absence of MSCs. The plot helps identify clusters and patterns of similarity or dissimilarity between samples.}
#'   \item{eigenvalues}{a barplot showing the percentage of explained variances by each principal component. This plot provides insights into the contribution of each principal component to the overall variation in the data.}
#'   \item{clustnames}{a  A list of cluster names with the highest contribution to PCA. This list helps identify the MSC clusters that have the most influence on the PCA results.}
#' @examples
#' data(matrices)
#' data(exData)
#' 
#' ### run function with all samples
#' res.pca <- lapply(matrices, function(x) msc.pca(x, samples = exData$samples, 
#'                   groups = exData$species, n=30, labels=FALSE, title=NULL))
#' 
#' res.pca$id95$eigenvalues
#' res.pca$id95$plot
#' 
#' ### use clusters with highest contribution to visualize in a heatmap
#' msc.heatmap(matrices[["id95"]][res.pca$id95$clustnames,], samples = exData$samples,
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

  
  #############   0   Tests   #############
  
  
  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (length(samples) != length(groups)) 
    stop("ERROR: samples and groups are not of equal length")
  if (length(samples) != ncol(clustmatrix)) {
    answer <- utils::menu(c("Yes", "No"), 
                   title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
    if (answer!=1) stop("Function stopped")
  }
  
  
  ### extra toevoegen: PC kiezen

  
  #############   1   PCA   #############
  
  
  pca <- FactoMineR::PCA(t(clustmatrix[,samples]),
                         scale.unit = F, 
                         ncp = 3, 
                         graph = F)

  
  
  if (labels == F) {
    plt <- factoextra::fviz_pca_ind(pca, 
                                    axes = c(1,2),
                                    geom.ind='point', 
                                    col.ind = 'black',
                                    fill.ind = groups,
                                    pointshape = 21,
                                    addEllipses = F, 
                                    legend.title = "", 
                                    alpha.ind = 0.8,
                                    pointsize = 5, 
                                    invisible = "quali", 
                                    title = title)
  } else {
    plt <- factoextra::fviz_pca_ind(pca, 
                                    axes = c(1,2),
                                    col.ind = 'black',
                                    fill.ind = groups,
                                    pointshape = 21,
                                    addEllipses = F, 
                                    legend.title = "", 
                                    alpha.ind = 0.8,
                                    pointsize = 5, 
                                    invisible = "quali", 
                                    title = title, 
                                    repel = TRUE,
                                    labelsize = 3)
  }

  contr <- factoextra::get_pca_var(pca)$contrib

  cl_index <- unique(order(contr[,1], decreasing = T)[1:n],
                     order(contr[,2], decreasing = T)[1:n],
                     order(contr[,3], decreasing = T)[1:n])
  cl_names <- rownames(clustmatrix[cl_index,samples])
  
  results <- list("plot" = plt,
                  "eigenvalues" = factoextra::fviz_eig(pca, addlabels = TRUE),
                  "clustnames" = cl_names)

  
  #############   2   Return PCA results    #############
  
  
  return(results)

  
  }

