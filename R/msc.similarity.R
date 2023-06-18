#' Minicircle Sequence Classes similarity
#'
#' The function msc.similarity returns a measure of minicircle sequence composition within and between groups of samples. Specifically, it estimates the absolute and relative number of Minicircle Sequence Classes (MSCs) that are unique to each group or shared between two or more groups. The function returns tables and barplots that summarize the number of unique or shared MSCs for each minimum percent identity (MPI) separately or combined over all MPIs.
#' 
#' @usage msc.similarity(clustmatrices, samples, groups)
#' @param clustmatrices a list of cluster matrices.
#' @param samples a vector containing the names of the samples. This can include all samples or it can be a subset. 
#' @param groups a vector, of equal length as samples,  specifying to which group (e.g. species) the samples belong to.
#' @return
#'   \item{absfreq}{a list per percent identity containing absolute frequency values of shared and unique MSCs.}
#'   \item{absfreq.plot}{a list of barplots visualizing previous results.}
#'   \item{relfreq}{a list per percent identity containing relative frequency values of shared and unique MSCs.}
#'   \item{relfreq.plot}{one barplot visualizing previous results.}
#' @examples
#' require(viridis)
#' data(matrices)
#' data(exData)
#' 
#' ### run function
#' sim <- msc.similarity(matrices, samples = exData$samples, 
#'                       groups = exData$species)
#' 
#' ### visualize results (absolute frequencies)
#' barplot(sim$absfreq$id93)
#' 
#' ### adjust plot (relative frequencies)
#' sim$relfreq.plot + scale_fill_viridis(discrete = TRUE)
#' 
#' sim$relfreq$id97["2"]*100 
#' sim$relfreq$id97["3"]*100 
#' 
#' ### reduce number of groups 
#' groups <- exData$species
#' levels(groups)[levels(groups)!='hybrid'] <- "non-hybrid"
#' sim.red <- msc.similarity(matrices, samples = exData$samples, groups = groups)
#' sim.red$relfreq.plot + scale_fill_viridis(discrete = TRUE)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom utils menu
#' @export

msc.similarity <- function(clustmatrices, samples, groups) {
  
  
  #############   0   Tests   #############
  
  
  if (!is.factor(groups)) stop("groups should be a factor")
  if (!is.list(clustmatrices)) stop("clustmatrices should be a list of matrices")
  if (length(samples) != length(groups)) stop("samples and groups are not of equal length")
  if (length(samples) != ncol(clustmatrices[[1]])) {
    answer <- utils::menu(c("Yes", "No"), 
                          title = "WARNING: You entered a subset of your samples.\nDo you wish to proceed?")
    if (answer != 1) stop("Function stopped")
  }
  
  
  #############   0     Define global functions or variables   #############
  
  
  Var1 <- Freq <- Var2 <- value <- NULL
  
  
  #############   1     Absolute frequency   #############
  
  
  id <- names(clustmatrices)
  num_levels <- length(levels(groups))
  
  absfreq <- absfreq.plots <- list()
  
  for (i in seq_along(id)) {
    temp <- list()
    for (s in seq_len(num_levels)) {
      ssp <- which(groups == levels(groups)[s])
      temp[[s]] <- rowSums(clustmatrices[[i]][, samples[ssp], drop = FALSE] != 0) / length(ssp)
    }
    tmp <- do.call("rbind", temp)
    absfreq[[i]] <- table(apply(tmp, 2, function(x) paste(which(x != 0), collapse = '')))
    names(absfreq[[i]])[names(absfreq[[i]]) == ""] <- 0
    
    df <- data.frame(absfreq[[i]])
    
    absfreq.plots[[i]] <- ggplot(df, aes(x = Var1, y = Freq)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(
        axis.text = element_text(hjust = 1),
        axis.title = element_text(face = "bold"),
        plot.caption = element_text(face = 'italic', color = "darkgray")
      ) +
      xlab('groups') + ylab('absolute MSC frequency') +
      geom_text(aes(label = Freq), vjust = -0.3, size = 3.5, color = "blue", fontface = "bold") +
      labs(caption = paste(seq_len(num_levels), ": ", levels(groups), "\n", collapse = '')) +
      ylim(0, max(absfreq[[i]]) + 100) +
      ggtitle(paste0("N of shared and unique MSC with ", id[i]))
  }
  
  
  #############   2     Relative frequency   #############
  
  
  relfreq <- lapply(absfreq, function(x) x / sum(x))
  grps <- unique(unlist(sapply(relfreq, names)))
  shared <- matrix(0, ncol = length(grps), nrow = length(relfreq), dimnames = list(id, grps))
  
  for (i in seq_len(length(relfreq))) {
    for (n in seq_len(length(grps))) {
      res <- relfreq[[i]][grps[n]]
      if (!is.null(res)) shared[i, grps[n]] <- res
    }
  }
  
  shared[is.na(shared)] <- 0
  
  relfreq.plot <- ggplot(reshape2::melt(shared), 
                         aes(x = as.factor(Var1), fill = as.factor(Var2), y = value)) +
    geom_bar(stat = "identity") +
    ylab('relative frequency') +
    ggtitle("N of shared and unique MSC across different % id's") +
    xlab("minimum percent identity") +
    theme_minimal() +
    labs(caption = paste(seq_len(num_levels), ": ", levels(groups), "\n", collapse = '')) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.caption = element_text(face = 'italic', color = "darkgray"),
      legend.title = element_blank()
    )
  
  names(relfreq) <- names(absfreq) <- names(absfreq.plots) <- id
  
  
  #############   3     Return   #############
  
  
  return(list(
    "absfreq" = absfreq,
    "absfreq.plots" = absfreq.plots,
    "relfreq" = relfreq,
    "relfreq.plot" = relfreq.plot
  ))
}