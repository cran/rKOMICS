#' Length of minicircles
#'
#' The msc.length function allows you to check the length of minicircle sequences based on a single FASTA file. This function helps determine the size distribution of minicircle sequences.
#' 
#' @usage msc.length(file, samples, groups)
#' @param file the name of the FASTA file that contains all the minicircle sequences. The file should be in the format "all.minicircles.circ.fasta".
#' @param samples a character vector containing the sample names.
#' @param groups a vector of the same length as the samples, specifying the groups (e.g., subspecies) to which the samples belong.
#' @return
#'   \item{length}{a numerical vector containing the lengths of the minicircle sequences. Each element corresponds to the length of a specific minicircle sequence.}
#'   \item{plot}{a histogram that visualizes the frequency distribution of minicircle sequence lengths. The histogram provides an overview of the length distribution of the minicircles.}
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


  #############   0   Tests   #############


  if (!file.exists(file)) stop("ERROR: File doesn't exist")
  if (length(file)<1) stop("ERROR: Your input parameter is empty")
  if (length(file)>1) stop("ERROR: Your input parameter is too long")


  #############   0   Define global functions or variables   #############
  
  
  freq.length <- NULL
  
  
  #############   1   Read in data    #############


  freq <- list()
  sequences <- ape::read.dna(file, 'fasta')
  samples <- sort(unique(gsub('_con.*','',attr(sequences, 'names'))))
  freq$length <- as.numeric(gsub('.*_len|_cir.*','',attr(sequences, 'names')))
  df <- data.frame(freq$length)
  
  freq$plot <- ggplot(df, aes(x=freq.length)) +
    geom_histogram(bins=50, alpha=0.8, show.legend = F, col="black", fill='gray') +
    xlab("Minicircle sequence length") + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold")) + 
    ylab("Frequency of minicircle length")


  #############   2   Return minicircle lengths    #############


  return(freq)

}
