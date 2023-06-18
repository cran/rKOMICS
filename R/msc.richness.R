#' Minicircle Sequence Cluster richness
#'
#' The msc.richness function calculates the measure of minicircle richness per sample by estimating the number of Minicircle Sequence Classes (MSCs) in each sample. It takes into account different minimum percent identities (MPIs) and returns a table of richness estimates per sample and per MPI. Additionally, it generates a boxplot that illustrates the minicircle richness across samples based on the estimated MSCs over a range of MPIs.
#'
#' @usage msc.richness(clustmatrices, samples, groups)
#' @param clustmatrices a list of cluster matrices obtained from the msc.matrix function. Each cluster matrix represents the presence or absence of MSCs in each sample.
#' @param samples a vector containing the names of the samples. 
#' @param groups a vector specifying the group (e.g., species) to which each sample belongs. 
#' @return
#'   \item{table}{a table summarizing the number of MSCs per sample at different percent identities. The table provides an overview of the estimated minicircle richness in each sample across the MPIs.}
#'   \item{plot}{a boxplot visualizing the minicircle richness across samples. The boxplot represents the distribution of richness estimates in each sample over the range of MPIs considered.}
#' @examples
#' require(ggplot2)
#' data(matrices)
#' data(exData)
#'
#' #### run function
#' richness <- msc.richness(matrices, samples = exData$samples, groups = exData$species)
#' 
#' apply(richness$table[which(richness$table$group=="L. peruviana"),-(1:2)], 2, mean)
#' apply(richness$table[which(richness$table$group=="L. braziliensis"),-(1:2)], 2, mean)
#' apply(richness$table[which(richness$table$group=="hybrid"),-(1:2)], 2, mean)  
#' 
#' #### visualize results
#' barplot(richness$table[,"id93"], names.arg = richness$table[,1],
#'         las=2, cex.names=0.4, main="N of MSC at id 93")
#' 
#' #### adjust plot
#' richness$plot + ggtitle("MSC richness across % id") + 
#'                 theme(axis.text.x = element_text(angle=45, hjust=1))
#' 
#' ### show results of subset 
#' table(exData$species)
#' hybrid <- which(exData$species=="hybrid")
#' # richness.subset <- msc.richness(matrices, samples = exData$samples[hybrid], 
#' #                                 groups = exData$species[hybrid])
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom utils menu
#'
#' @export

msc.richness <- function(clustmatrices, samples, groups) {
  
  
  #############   0   Tests   #############
  
  
  if (!is.factor(groups)) stop("groups should be a factor")
  if (!is.list(clustmatrices)) stop("clustmatrices should be a list of matrices")
  if (length(samples) != length(groups)) stop("samples and groups are not of equal length")
  if (length(samples) != ncol(clustmatrices[[1]])) {
        answer <- utils::menu(c("Yes", "No"), 
                        title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
        if (answer!=1) stop("Function stopped")
  }
  
  
  #############   0     Define global functions or variables   #############
  
  
  sample <- group <- value <- NULL
  
  
  #############   1   MSC per sample across different % id    #############
  
  
  id <- names(clustmatrices)
  
  N_between <- data.frame(matrix(nrow = length(samples), ncol = length(id)+2))
  N_between[,1:2] <- c(samples,as.character(groups))
  colnames(N_between) <- c("sample", "group",id)
  
  for (i in 1:length(id)) {
    N_between[, i + 2] <- colSums(clustmatrices[[i]][, samples] > 0)
  }
  melted <- reshape2::melt(N_between, id=c("sample", "group"))
  
  plt <- ggplot(melted, aes(x=sample, y=value, col=group, fill=group)) + 
    geom_boxplot(alpha=0.5, show.legend = F) +
    xlab('') + 
    theme_minimal() +
    facet_grid(. ~ group, scales='free', space='free') + 
    ylab('Number of MSC across different % id') +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.title = element_text(face = "bold"))

  
  #############   2   Return MSC richness    #############
  
  
  return(list("table" = N_between, "plot" = plt))

  
}