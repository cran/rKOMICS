#' Check the read depth of assembled minicircles
#'
#' The depth statistics, generated with KOMICS, include average, median, minimum and maximum per site read depth of every minicircle contig that has been assembled.
#' The msc.depth function allows you to summarize those statistics and to estimate minicircle copy numbers by standardizing median read depths per minicircle contig to the median genome-wide read depths.
#'
#' @usage msc.depth(depthstats, groups, HCN = NULL)
#' @param depthstats a character vector containing the file names of the depth statistics (output of KOMICS), e.g. sampleA.depthstats.txt, sampleB.depthstats.txt,... .
#' @param groups a vector specifying to which groups (e.g. species) the samples belong to.
#' @param HCN a numeric vector containing haploid copy numbers of the corresponding samples (optional, by default set to null).
#' @return 
#'   \item{all}{a table merging depth statistics of all samples. Depth statistics include the average, median, minimum and maximum per site read depth. }
#'   \item{plots}{a plot per sample visualizing the median read depth distribution.}
#'   \item{medianRD}{one graph summarizing the median read depth distribution of all samples.}
#'   \item{CN}{one graph summarizing the copy number (if HCN is not null) of all samples.}
#' @examples
#' require(ggpubr)
#' data(exData)
#'
#' ### run function
#' \donttest{
#' depth <- msc.depth(depthstats = system.file("extdata", 
#'                   exData$depthstats, package = "rKOMICS"), groups = exData$species,
#'                   HCN = exData$medGWD/2)
#' 
#' ### visualize results 
#' hist(depth$all[,"MEDIAN.DEPTH"], breaks=100, 
#'     main="Global median depth distribution",xlab = (''))
#' 
#' ### alter plot 
#' annotate_figure(depth$plots$CUM29A1, fig.lab = "CUM29A1", 
#'                fig.lab.pos = "bottom.right", fig.lab.face = 'italic')
#' }
#' @importFrom utils read.table read.csv
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr group_by summarise
#' @import ggplot2
#' @export

msc.depth <- function(depthstats, groups, HCN = NULL) {
  
  
  ############# tests
  
  
  if (sum(file.exists(depthstats))!=length(depthstats)) stop("ERROR: One or more files doesn't exist")
  if (!is.null(HCN)) {
    if (length(depthstats) != length(HCN)) stop("ERROR: the HCN vector and depthstats vector should be of equal length.")
  }
  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (length(depthstats)!=length(groups)) stop('ERROR: groups factor and depthstats are not of equal length. ')
  
  
  ############# depth stats per minicircle contig
  
  
  theme_set(theme_minimal() + theme(axis.title = element_text(face = 'bold')))
  
  depths <- depth_plts <- hists_MEDIAN.DEPTH <- box_MEDIAN.DEPTH <- box_CN <- list()
  
  for (n in 1:length(depthstats)) {
    
    depths[[n]] <- utils::read.csv(depthstats[n], sep=' ', header=T, row.names = 1)[,1:4]
    depths[[n]][,"sample"] <- gsub("_contig.*","",rownames(depths[[n]])[1])
    depths[[n]][,"group"] <- groups[n]
    
    hists_MEDIAN.DEPTH[[n]] <- ggplot(depths[[n]], aes(x=depths[[n]]$MEDIAN.DEPTH)) +
                                   geom_histogram(bins = 50, color="black", fill="gray") + 
                                   ylab('frequency') + xlab('median depth') + ggtitle(depths[[n]][1,5])
    
    box_MEDIAN.DEPTH[[n]] <- ggplot(depths[[n]], aes(x=depths[[n]]$MEDIAN.DEPTH)) +
                                    geom_boxplot() + xlab('median depth') + 
                                    theme(axis.text.y = element_blank()) 

    if (!is.null(HCN)) {
      
      depths[[n]]$HCN <- HCN[n]
      depths[[n]]$CN <- depths[[n]]$MEDIAN.DEPTH/HCN[n]
      depths[[n]]$MIN.CN <- depths[[n]]$MIN.DEPTH/HCN[n]
      depths[[n]]$MAX.CN <- depths[[n]]$MAX.DEPTH/HCN[n]
      
      box_CN[[n]] <- ggplot(depths[[n]], aes(x=depths[[n]]$CN)) +
                                geom_boxplot() + xlab('copy number') + 
                                theme(axis.text.y = element_blank()) 

      depth_plts[[n]] <- ggpubr::ggarrange(hists_MEDIAN.DEPTH[[n]], box_CN[[n]], nrow=2, heights = c(3,1), align = "v")
  
    } else {
      
      depths[[n]]$HCN <- depths[[n]]$CN <- depths[[n]]$MIN.CN <- depths[[n]]$MAX.CN <- NA
      
      depth_plts[[n]] <- ggpubr::ggarrange(hists_MEDIAN.DEPTH[[n]], box_MEDIAN.DEPTH[[n]], nrow=2, heights = c(3,1), align = "v")
    }


  }
  
  depths <- do.call("rbind", depths)
  names(depth_plts) <- unique(depths$sample)
  
  mean.RD <- data.frame(depths %>% dplyr::group_by(sample) %>% dplyr::summarise(mean(depths$MEDIAN.DEPTH), .groups = 'drop'))
  
  box_MEDIAN.DEPTH <- ggplot(depths, aes(y=depths$MEDIAN.DEPTH, x=sample, fill=depths$group, color=depths$group)) +
    geom_boxplot(alpha=0.5, show.legend = F) + xlab('') + 
    theme(axis.text.x = element_text(angle=90, hjust=1)) + ylab('Minicircle median read depth') + 
    facet_grid(. ~ group, scales="free", space="free")
  
  if (!is.null(HCN)) {
    box_CN <- ggplot(depths, aes(y=depths$CN, x=sample, fill=depths$group, color=depths$group)) +
      geom_boxplot(alpha=0.5, show.legend = F) + xlab('') + 
      theme(axis.text.x = element_text(angle=90, hjust=1)) + ylab('Minicircle copy number') + 
      facet_grid(. ~ group, scales="free", space="free")
  } else {box_CN <- NA}
  
  
  ############# return
  
  
  return(list("all" = depths,
              "plots" = depth_plts,
              "medianRD" = box_MEDIAN.DEPTH,
              "CN" = box_CN))
  
}



