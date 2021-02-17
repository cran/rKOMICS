#' Retrieve sequences 
#'
#' The msc.seqs retrieves the DNA sequence of a minicircle sequence cluster (MSC) together with all its hit sequences. 
#'
#' @usage msc.seqs(fastafile, ucfile, clustnumbers, writeDNA = TRUE)
#' @param fastafile the name of the fasta file containing all minicircle sequences. 
#' @param ucfile the name of the uc file. 
#' @param clustnumbers a character vector containing the cluster numbers (in the format "C0", "C1", ...) of which cluster and hit sequences need to be extracted. 
#' @param writeDNA a logical parameter which is by default set to TRUE. This will write fasta files to the current directory. 
#' @return a table which summarizes the number of hit sequences found in the MSC, the MSC name and where the MSC is present (strain names). 
#' @return one fasta file per MSC with all its hits sequences. 
#' @examples
#' data(exData)
#'
#' ### select a subset of MSC
#' Lpe <- which(exData$species == "L. peruviana")
#' specific <- msc.subset(matrices[[7]], subset = Lpe)
#' 
#' ### run function
#' seq <- msc.seqs(fastafile = system.file("extdata", "all.minicircles.circ.fasta", package="rKOMICS"),
#'                 ucfile = system.file("extdata", exData$ucs, package="rKOMICS")[7], 
#'                 clustnumbers = specific$clustnumbers, writeDNA = FALSE)
#'
#' @importFrom ape write.dna
#' @importFrom stats na.omit
#' @export

msc.seqs <- function(fastafile, ucfile, clustnumbers, writeDNA = TRUE) {

  
  ############# tests
  
  
  if (!file.exists(fastafile)) stop("ERROR: fasta file doesn't exist")
  if (!file.exists(ucfile)) stop("ERROR: uc file doesn't exist")
  if (length(clustnumbers) == 0) stop("ERROR: clustnumbers vector is empty")
  if (!is.vector(clustnumbers)) stop("ERROR: clustnumbers should be a vector ")
  if (!is.logical(writeDNA)) stop("ERROR: writeDNA should be logical")
  
  ############# read in all sequences
  
  
  sequences <- ape::read.dna(fastafile, 'fasta')
  
  
  ############# read in uc
  
  
  id <- as.numeric(gsub(".uc", "", gsub(".*id", "", ucfile)))

  uc_H <- read.uc(ucfile)$hits
  uc_C <- read.uc(ucfile)$clusters
  uc_C2 <- read.uc(ucfile)$clustnumbers

  
  ############# retrieve sequences
  
  
  summary <- data.frame(matrix(nrow=length(clustnumbers), ncol=3))
  rownames(summary) <- clustnumbers
  colnames(summary) <- c("N contigs", "cluster strain", "hit strains")

  for (n in 1:length(clustnumbers)) {
    dat <- uc_H[which(paste('C', uc_H$V2, sep='') == clustnumbers[n]),]
    datcontigs <- unique(c(as.character(dat$V9), as.character(dat$V10)))
    dna2 <- sequences[datcontigs]
    if (length(dna2) != 0) {
      #dna2 <- dna2[-which(sapply(dna2, is.null))]
      names(dna2) <- gsub('_contig.*','',names(dna2))
      if (writeDNA == TRUE) {
          ape::write.dna(x = dna2, file = paste(clustnumbers[n],'.fasta',sep=''), format = 'fasta')}
      summary[n, ] <- c(length(datcontigs),
                        dat$V10[1],
                        paste(gsub('_contig.*','',datcontigs), sep="", collapse=", "))
    }

  }

  return(list("summary" = stats::na.omit(summary)))

}
