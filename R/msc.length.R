#' Length of minicircles
#'
#' The msc.length function allows you to check the length of minicircle sequences based on one fasta file. 
#' 
#' @usage msc.length(file, samples, groups)
#' @param file the name of the fasta file containing all minicircle sequences, e.g. all.minicircles.circ.fasta.
#' @param samples a character vector containing the sample names.
#' @param groups a vector, of equal length as samples, specifying to which group (e.g. subspecies) the samples belong to.
#' @return
#'   \item{length}{a numerical vector containing the length of minicircle sequences.}
#'   \item{plot}{a histogram visualizing the length frequency of minicircle sequences.}
#' @examples
#' require(ggplot2)
#' require(ggpubr)
#' 
#' ### run function
#' bf <- msc.length(file = system.file("extdata", "all.minicircles.fasta", package="rKOMICS"),
#'                  samples = exData$samples, groups = exData$subspecies)
#' af <- msc.length(file = system.file("extdata", "all.minicircles.circ.fasta", package="rKOMICS"),
#'                  samples = exData$samples, groups = exData$subspecies)
#' 
#' length(which(bf$length<800)) 
#' length(which(bf$length>1400)) 
#' 
#' ### visualize results
#' hist(af$length, breaks=50)
#' 
#' ### alter plot
#' ggarrange(bf$plot + labs(caption = "Before filtering"), 
#'           af$plot + labs(caption = "After filtering"), nrow=2)
#' 
#' 
#' @import ggplot2
#' @importFrom ape read.dna
#' @export

msc.length <- function(file, samples, groups) {


  ############# tests


  if (!file.exists(file)) stop("ERROR: File doesn't exist")
  if (length(file)<1) stop("ERROR: Your input parameter is empty")
  if (length(file)>1) stop("ERROR: Your input parameter is too long")


  ############# read in data


  freq <- list()
  sequences <- ape::read.dna(file, 'fasta')
  samples <- sort(unique(gsub('_con.*','',attr(sequences, 'names'))))
  freq$length <- as.numeric(gsub('.*_len|_cir.*','',attr(sequences, 'names')))
  df <- data.frame(freq$length)
  freq$plot <- ggplot(df, aes(x=df$freq.length)) +
    geom_histogram(bins=50, alpha=0.8, show.legend = F, col="black", fill='gray') +
    xlab("Minicircle sequence length") + theme_minimal() +
    theme(axis.title = element_text(face = "bold")) + ylab("Frequency of minicircle length")


  ############# return


  return(freq)

}
