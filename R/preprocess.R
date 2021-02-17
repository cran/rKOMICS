#' Filtering of minicircle sequences
#'
#' Assembling minicircle sequences with KOMICS generates individual fasta files (one per sample). 
#' The preprocess function allows you to filter the minicircle sequences based on sequence length (as the size of minicircular kDNA is species-specific and variable) and circularization success. 
#' The function will write filtered individual fasta files in the current working directory.
#'
#' @usage preprocess(files, groups, circ = TRUE, min = 500, max = 1500, writeDNA = TRUE)
#' @param files a character vector containing the fasta file names in the format sampleA.minicircles.fasta, sampleB.minicircles.fasta,... (output of KOMICS).
#' @param groups a factor specifying to which group (e.g. species) the samples belong to. It should have the same length as the list of files.
#' @param circ a logical parameter. By default non-circularized minicicle sequences will be excluded. If interested in non-circularized sequences as well, set the parameter to FALSE.
#' @param min a minimum value for the minicircle sequences length. Default value is set to 500.
#' @param max a maximum value for the minicircle sequences length. Default value is set to 1500.
#' @param writeDNA a logical parameter. By default filtered minicircle sequences will by written in fasta format to the current working directory. Set to FALSE if only interested in other output values like plots and summary.
#' @return
#'   \item{samples}{the sample names (based on the input files).}
#'   \item{N_MC}{a table containing the sample name, which group it belongs to and the number of minicirce sequences (N_MC) before and after filtering.}
#'   \item{plot}{a barplot visualizing the number of minicircle sequences per sample before and after filtering.}
#'   \item{summary}{the total number of minicircle sequences before and after filtering.}
#' @examples
#' require(ggplot2)
#' data(exData)
#'
#' ### setwd("")
#' 
#' ### run function
#' table(exData$species)
#' pre <- preprocess(files = system.file("extdata", exData$fastafiles, package="rKOMICS"),
#'                   groups = exData$species,
#'                   circ = TRUE, min = 500, max = 1200, writeDNA = FALSE)
#' 
#' pre$summary 
#' 
#' ### visualize results
#' barplot(pre$N_MC[,"beforefiltering"], 
#'         names.arg = pre$N_MC[,1], las=2, cex.names=0.4)
#' 
#' ### alter plot
#' pre$plot + labs(caption = paste0('N of MC sequences before and after filtering, ', Sys.Date()))
#'
#' @import ggplot2
#' @importFrom ape read.dna
#' @export

preprocess <- function(files, groups, circ = TRUE, min = 500, max = 1500, writeDNA = TRUE) {


  ############# tests


  if (sum(file.exists(files))!=length(files)) stop("ERROR: One or more files doesn't exist")
  if (!is.factor(groups)) stop("ERROR: groups should be a factor")
  if (length(groups) != length(files)) stop("ERROR: Input files and groups should be of equal length")
  if (!is.logical(c(circ, writeDNA))) stop("ERROR: circ/writeDNA should be logical")
  if (!is.numeric(c(min, max)) & length(c(min, max))!=2) stop("ERROR: min/max should be one numeric value")


  ############# filter MC-sequences


  samples <- sort(gsub(".*/", "", gsub(".minicircles.fasta", "", files)))
  pp <- data.frame(samples=samples, groups=groups)

  for (i in 1:length(files)) {
    dna <- ape::read.dna(files[i], 'fasta')
    pp$beforefiltering[i] <- length(dna)
    dna <- dna[which(as.numeric(gsub('.*_len|_cir.*','',attr(dna, 'names')))<max &
                       as.numeric(gsub('.*_len|_cir.*','',attr(dna, 'names')))>min)]
    if (circ == TRUE) {
      dna <- dna[grep("circ", attr(dna, 'names'))]
    }
    
    if (length(dna)!=0 && writeDNA == TRUE) {
      ape::write.dna(dna, paste0(pp$sample[i], ".minicircles.circ.fasta"), "fasta")
    }
    
    pp$afterfiltering[i] <- length(dna)
  }


  ############# plot

  
  gplot <- ggplot(pp, aes(y=pp$beforefiltering, x=samples)) +
    geom_bar(stat="identity", alpha=0.3, color="gray90", show.legend = F) + ylab('Number of minicircle sequences') +
    facet_grid(. ~ groups, scales='free', space='free') + theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.title = element_text(face="bold")) + xlab("") +
    geom_bar(aes(y=pp$afterfiltering, x=samples, fill=groups), stat="identity", alpha=0.7, show.legend = F)


  ############# return


  results <- list(samples, pp, gplot, colSums(pp[,3:4]))
  names(results) <- c("samples", "N_MC", "plot", "summary")
  return(results)

}
