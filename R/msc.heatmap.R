#' Visualization of cluster matrices
#'
#' The msc.heatmap function generates a heatmap that summarizes the presence or absence of Minicircle Sequence Classes (MCSs) between groups of samples. It takes an input cluster matrix, generated using the msc.matrix function, and visualizes the clustering patterns of MCSs.
#'
#' @usage msc.heatmap(clustmatrix, samples, groups)
#' @param clustmatrix a cluster matrix generated with the msc.matrix function. This matrix represents the presence or absence of MCSs in each sample.
#' @param samples a vector containing the sample names. 
#' @param groups a vector specifying the groups (e.g., species) to which the samples belong.
#' @return a heatmap that visualizes the clustering patterns of MCSs between sample groups. The heatmap provides an overview of the presence or absence of MCSs and helps identify shared or distinct MCS patterns among the groups.
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


  #############   0   Tests   #############


  if (!is.matrix(clustmatrix)) stop("clustmatrix should be a matrix")
  if (length(samples) != ncol(clustmatrix)) {
    answer <- utils::menu(c("Yes", "No"), 
                          title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
    if (answer!=1) stop("Function stopped")
  }
  
  
  #############   1   Heatmap   #############

  
  groups <- as.factor(groups)
  ncols <- length(unique(groups))+1
  col_fun <- circlize::colorRamp2(c(0, 1, 2, max(clustmatrix[,samples])), 
                                  c("white","pink", "red", "black"))
  
  captions <- paste0(seq_along(levels(groups)), ": ", levels(groups))
  
  
  #############   2   Calculate MSCgroups   #############
  
  
  temp <- lapply(seq_along(levels(groups)), function(s) {
    ssp <- which(groups == levels(groups)[s])
    apply(clustmatrix[, ssp, drop = FALSE], 1, function(x) sum(x != 0)) / length(ssp)
  })
  
  tmp <- do.call("rbind", temp)
  MSCgroups <- apply(tmp, 2, function(x) paste(which(x != 0), collapse=''))
  
  
  #############   3   Create Heatmap    #############
  
  
  heatmap <- ComplexHeatmap::Heatmap(t(clustmatrix[, samples]),
                                     col = col_fun,
                                     show_column_names = FALSE,
                                     column_names_gp = grid::gpar(fontsize = 2),
                                     column_title_gp = grid::gpar(fontsize = 7),
                                     column_title_rot = 90,
                                     column_split = MSCgroups,
                                     column_km = length(levels(MSCgroups)),
                                     row_title = "samples",
                                     row_names_gp = grid::gpar(fontsize = 7),
                                     name = "MSC",
                                     show_column_dend = TRUE)
  
  
  #############   4   Add Groups Heatmap if applicable    #############
  
  
  if (length(unique(groups)) != 1) {
    groups_heatmap <- ComplexHeatmap::Heatmap(groups,
                                              name = " ",
                                              width = grid::unit(5, "mm"),
                                              col = 2:ncols)
    heatmap <- ComplexHeatmap::draw(heatmap + groups_heatmap,
                                    auto_adjust = FALSE,
                                    column_title = paste(captions, collapse = "  -  "),
                                    column_title_side = "bottom",
                                    column_title_gp = grid::gpar(fontsize = 10))
  } else {
    heatmap <- ComplexHeatmap::draw(heatmap,
                                    auto_adjust = FALSE,
                                    column_title = paste(captions, collapse = "  -  "),
                                    column_title_side = "bottom",
                                    column_title_gp = grid::gpar(fontsize = 10))
  }
  
  
  #############   5   Return the Heatmap    #############
 
  
  return(heatmap)

  
  }