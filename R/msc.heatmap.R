#' Visualization of cluster matrices
#'
#' The msc.heatmap function generates a heatmap of the input cluster matrix that summarizes the presence or absence  of Minicircle Cluster Sequences (MCSs) between groups of samples.
#'
#' @usage msc.heatmap(clustmatrix, samples, groups)
#' @param clustmatrix one cluster matrix generated with the msc.matrix function.
#' @param samples a vector containing the sample names. This can include all samples or it can be a subset. 
#' @param groups a vector specifying to which groups (e.g. species) the samples belong to.
#' @return a heatmap
#' @examples
#' data(exData)
#' data(matrices)
#'
#' ### run function 
#' msc.heatmap(matrices[["id80"]], groups = exData$species,
#'             samples = exData$samples )
#'    
#' ### run function on every cluster matrix with subset of samples
#' ### you will be asked to confirm
#' table(exData$species)
#' hybrid <- which(exData$species=="hybrid")
#' # msc.heatmap(matrices[["id97"]], groups = exData$species[hybrid], 
#' #             samples = exData$samples[hybrid])
#' 
#' 
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid gpar
#' @importFrom utils menu
#' @export

msc.heatmap <- function(clustmatrix, samples, groups) {


  ############# tests


  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (!is.matrix(clustmatrix)) stop("ERROR: clustmatrix should be a matrix")
  if (length(samples) != ncol(clustmatrix)) {
    answer <- utils::menu(c("Yes", "No"), 
                   title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
    if (answer!=1) stop("Function stopped")
  }
  
  
  ############# msc.heatmap

  
  ncols <- length(unique(groups))+1
  col_fun = circlize::colorRamp2(c(0,1,2,max(clustmatrix[,samples])), c("white","pink", "red", "black"))
  
  captions <- vector()
  for (i in 1:length(levels(groups))) {
    captions[i] <- paste0(i, ": ", levels(groups)[i])
  }
  
  
  temp <- list()
  for (s in 1:length(levels(groups))) {
    ssp <- which(groups==levels(groups)[s])
    temp[[s]] <- apply(clustmatrix[,ssp, drop=F], 1, function(x) sum(x!=0))/length(ssp)
  }
  tmp <- do.call("rbind", temp)
  MSCgroups <- apply(tmp, 2, function(x) paste(which(x != 0), collapse=''))
  
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(t(clustmatrix[,samples]), col=col_fun, show_column_names = F,
                                               column_names_gp = grid::gpar(fontsize = 2),
                                               column_title_gp = grid::gpar(fontsize = 7),
                                               column_split = MSCgroups,
                                               column_km = length(levels(MSCgroups)),
                                               row_title = "samples", 
                                               row_names_gp = grid::gpar(fontsize = 7), name="MSC",
                                               show_column_dend = T ) +
                         if (length(unique(groups)) != 1) {
                           ComplexHeatmap::Heatmap(groups, name = " ", width = grid::unit(5, "mm"), 
                                                   col=c(2:ncols))
                         },
                       auto_adjust=F, column_title = paste(captions,collapse="  -  "),
                       column_title_side = "bottom", column_title_gp = grid::gpar(fontsize = 10)) 
  
}


  