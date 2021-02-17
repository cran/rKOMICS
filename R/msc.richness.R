#' Minicircle Sequence Cluster richness
#'
#' The msc.richness function counts how many Minicircle Sequence Clusters (MSC) are present per sample across different percent identities.
#'
#' @usage msc.richness(clustmatrices, samples, groups)
#' @param clustmatrices a list of cluster matrices.
#' @param samples a vector containing the names of the samples. This can include all samples or it can be a subset. 
#' @param groups a vector, of equal length as samples, specifying to which group (e.g. species) the samples belong to.
#' @return
#'   \item{table}{a table containing the number of MSC per sample across different percent identities.}
#'   \item{plot}{a boxplot visualizing previous results. }
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
  
  
  ############# tests
  
  
  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (!is.list(clustmatrices)) stop("ERROR: clustmatrices should be a list of matrices")
  if (length(samples) != length(groups)) stop("ERROR: samples and groups are not of equal length")
  if (length(samples) != ncol(clustmatrices[[1]])) {
        answer <- utils::menu(c("Yes", "No"), 
                        title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
        if (answer!=1) stop("Function stopped")
  }
  
  
  ############# MSC per sample across different % id
  
  
  id <- names(clustmatrices)
  
  N_between <- data.frame(matrix(nrow = length(samples), ncol = length(id)+2))
  N_between[,1:2] <- c(samples,as.character(groups))
  colnames(N_between) <- c("samples", "group",id)
  
  for (i in 1:length(id)) {
    N_between[,i+2] <- apply(clustmatrices[[i]][,samples], 2, function(x) sum(x>0))
  }
  melted <- reshape2::melt(N_between, id=c("samples", "group"))
  
  plt <- ggplot(melted, aes(x=samples, y=melted$value, col=melted$group, fill=melted$group)) + geom_boxplot(alpha=0.5, show.legend = F) +
    xlab('') + theme_minimal() +
    facet_grid(. ~ group, scales='free', space='free') + ylab('Number of MSC across different % id') +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.title = element_text(face = "bold"))

  
  ############# return
  
  
  return(list("table" = N_between, "plot" = plt))

  
}