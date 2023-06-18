#' Filtering of minicircle sequences
#'
#' The preprocess function is used to filter minicircle sequences based on sequence length and circularization success. When minicircle sequences are assembled with KOMICS, individual fasta files are generated for each sample. This function allows you to filter these sequences based on their length and whether they are circularized or not. The filtered sequences are then written into individual FASTA files in the current working directory.
#'
#' @usage preprocess(files, groups, circ = TRUE, min = 500, max = 1500, writeDNA = TRUE)
#' @param files a character vector containing the names of the fasta files. Each file corresponds to the minicircle sequences of a specific sample. The file names should be in the format sampleA.minicircles.fasta, sampleB.minicircles.fasta, and so on (output of KOMICS).
#' @param groups a factor specifying the group (e.g., species) to which each sample belongs. It should have the same length as the list of files, indicating the group assignment for each sample.
#' @param circ a logical parameter that determines whether non-circularized minicircle sequences should be included or excluded from the filtering process. By default, non-circularized sequences are excluded (circ = TRUE). If you are interested in including non-circularized sequences, you can set the parameter to FALSE.
#' @param min the minimum length threshold for filtering minicircle sequences. Sequences with a length below this threshold will be excluded. The default value is set to 500.
#' @param max the maximum length threshold for filtering minicircle sequences. Sequences with a length above this threshold will be excluded. The default value is set to 1500.
#' @param writeDNA a logical parameter that determines whether the filtered minicircle sequences should be written in FASTA format to the current working directory. By default, the filtered sequences are written (writeDNA = TRUE). If you are only interested in other output values like plots and summary, you can set this parameter to FALSE.
#' @return
#'   \item{samples}{the sample names based on the input files.}
#'   \item{N_MC}{a table containing the sample names, the corresponding group assignment, and the number of minicircle sequences (N_MC) before and after filtering.}
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


  #############   0   Tests   #############


  if (sum(file.exists(files))!=length(files)) stop("One or more files doesn't exist")
  if (length(groups) != length(files)) stop("Input files and groups should be of equal length")
  if (!is.logical(c(circ, writeDNA))) stop("circ/writeDNA should be logical")
  if (!is.numeric(c(min, max)) & length(c(min, max))!=2) stop("min/max should be one numeric value")

  
  #############   0   Define global functions or variables   #############
 
  
   beforefiltering <- afterfiltering <- NULL
  
  
  #############   1   Filter MC-sequences   #############


  samples <- sort(gsub(".*/|\\.minicircles.fasta", "", files))
  pp <- data.frame(samples=samples, groups=groups)
  groups <- as.factor(groups)
  
  for (i in seq_along(files)) {
    dna <- ape::read.dna(files[i], 'fasta')
    pp$beforefiltering[i] <- length(dna)
    dna <- dna[as.numeric(gsub('.*_len|_cir.*', '', attr(dna, 'names'))) < max &
                 as.numeric(gsub('.*_len|_cir.*', '', attr(dna, 'names'))) > min]
    if (circ) {
      dna <- dna[grep("circ", attr(dna, 'names'))]
    }
    
    if (length(dna) != 0 && writeDNA) {
      ape::write.dna(dna, paste0(pp$samples[i], ".minicircles.circ.fasta"), "fasta")
    }
    
    pp$afterfiltering[i] <- length(dna)
  }


  #############   2   Plot   #############


  gplot <- ggplot(pp, aes(y = beforefiltering, x = samples)) +
    geom_bar(stat = "identity", alpha = 0.3, color = "gray90", show.legend = FALSE) +
    ylab("Number of minicircle sequences") +
    facet_grid(. ~ groups, scales = 'free', space = 'free') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title = element_text(face = "bold")) +
    xlab("") +
    geom_bar(aes(y = afterfiltering, x = samples, fill = groups), 
             stat = "identity", alpha = 0.7, show.legend = FALSE)

  
  #############   3   Return   #############


  results <- list(samples, pp, gplot, colSums(pp[,3:4]))
  names(results) <- c("samples", "N_MC", "plot", "summary")
  return(results)

}
