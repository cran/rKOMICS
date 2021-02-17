#' Minicircle Sequence Cluster similarity
#'
#' The msc.similarity function allows you to check the absolute and relative frequency of shared and unique MSC between different groups across different percent identities.
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


  ############# tests


  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (!is.list(clustmatrices)) stop("ERROR: clustmatrices should be a list of matrices")
  if (length(samples) != length(groups)) stop("ERROR: samples and groups are not of equal length")
  if (length(samples) != ncol(clustmatrices[[1]])) {
    answer <- utils::menu(c("Yes", "No"), 
                   title="WARNING: You entered a subset of your samples.\nDo you wish to procede?")
    if (answer!=1) stop("Function stopped")
  }


  ############# abs frequency

  
  id <- names(clustmatrices)


  N_shared <- matrix(ncol = length(samples), nrow = length(id))
  colnames(N_shared) <- samples
  relfreq <- absfreq <- temp <- list()

  for (i in 1:length(id)) {
    for (s in 1:length(levels(groups))) {
      ssp <- which(groups==levels(groups)[s])
      temp[[s]] <- apply(clustmatrices[[i]][,ssp, drop=F], 1, function(x) sum(x!=0))/length(ssp)
    }
    tmp <- do.call("rbind", temp)
    absfreq[[i]] <- table(apply(tmp, 2, function(x) paste(which(x != 0), collapse='')))
  }

  captions <- vector()
  for (i in 1:length(levels(groups))) {
    captions[i] <- paste0(i, ": ", levels(groups)[i], "\n")
  }

  absfreq.plots <- list()

  for (i in 1:length(absfreq)) {
    df <- data.frame(absfreq[[i]])
    absfreq.plots[[i]] <- ggplot(df, aes(x=df$Var1, y=df$Freq)) +
      geom_bar(stat = "identity") + theme_minimal() +
      theme(axis.text = element_text(hjust = 1),
            axis.title = element_text(face="bold"),
            plot.caption = element_text(face='italic', color="darkgray")) +
      xlab('groups') + ylab('absolute MSC frequency') +
      geom_text(aes(label=df$Freq), vjust=-0.3, size=3.5, color="blue", fontface="bold") +
      labs(caption = paste(captions, collapse='')) + ylim(0,max(absfreq[[i]])+100) +
      ggtitle(paste0("N of shared and unique MSC with a % id of ", id[i]))
  }


  ############# rel frequency


  relfreq <- lapply(absfreq, function(x) x/sum(x))
  grps <- names(relfreq[[which.max(sapply(relfreq, length))]])


  shared <- matrix(ncol = length(grps),nrow=length(relfreq))
  colnames(shared) <- grps
  rownames(shared) <- id
  for (i in 1:length(relfreq)) {
    for (n in 1:length(grps)) {
      res <- relfreq[[i]][which(names(relfreq[[i]]) == grps[n])]
      if (length(res) > 0) shared[i,grps[n]] <- res
    }
  }

  ### test:
  #apply(shared, 1, function(x) sum(x, na.rm = T))

  shared[is.na(shared)] <- 0

  melted <- reshape2::melt(shared)
  colnames(melted) <- c("% id","group", "% shared")
  melted$group <- as.factor(melted$group)
  melted$`% id` <- as.factor(melted$`% id`)
  relfreq.plot <- ggplot(melted, aes(x=melted$`% id`,fill=melted$group, y=melted$`% shared`)) +
    geom_bar(stat="identity") +
    ylab('relative frequency') +
    ggtitle("N of shared and unique MSC across different % id's") +
    theme_minimal() +
    labs(caption = paste(captions, collapse='')) + 
    theme(axis.title = element_text(face="bold"),
          plot.caption = element_text(face='italic', color="darkgray"))


  names(relfreq) <- names(absfreq) <- names(absfreq.plots) <- id
  
  
  ############# return


  return(list("absfreq" = absfreq, "absfreq.plots" = absfreq.plots,
              "relfreq" = relfreq, "relfreq.plot" = relfreq.plot))

}
