#' Cluster Analyses
#'
#' The msc.uc function allows you to check some clustering : insertions and deletions --> ask Fre.
#' The total amount of minicircle sequence clusters (MSCs) found per percent identity will be calculated as well as the number of gaps.
#'
#' @usage msc.uc(files)
#' @param files a character vector containing the uc file names (output of VSEARCH) e.g. all.minicircles.circ.id70.uc, all.minicircles.circ.id80.uc...
#' @return
#'   \item{MSCs}{a numerical vector containing the number of MSC per percent identity.}
#'   \item{perfect aligments}{a numerical vector containing the proportions of perfect alignments per percent identity.}
#'   \item{insertions}{a table showing the length and the number of insertions across different percent identities.}
#'   \item{deletions}{a table showing the length and the number of deletions across different percent identities.}
#'   \item{plots}{different plots showing previous results.}
#' @examples
#' data(exData)
#' 
#' ### run function
#' \donttest{
#' ucs <- msc.uc(files = system.file("extdata", exData$ucs, package="rKOMICS"))
#' 
#' ucs$MSCs["100"] 
#' ucs$MSCs["97"] 
#' 
#' ### results
#' ucs$plots
#' }
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stringr str_count
#' @importFrom ggpubr ggarrange annotate_figure
#' @export


msc.uc <- function(files) {


  ############# tests


  if (sum(file.exists(files))!=length(files)) stop("One or more files doesn't exist")
  
  
  ############# analysis 
  
  
  id <- as.numeric(regmatches( gsub(".*/", "", files), gregexpr("[[:digit:]]+",  gsub(".*/", "", files))))
  
  
  MSCs <- PROP.PA <- vector()
  gplots <- pltINS <- pltDEL <- list()
  INS <- DEL <- matrix(ncol = 8, nrow = length(id))
  colnames(INS) <- c('1bp','2bp','3bp','minimum 4bp','1 insertion','2 insertions','3 insertions', 'minimum 4 insertions')
  colnames(DEL) <- c('1bp','2bp','3bp','minimum 4bp','1 deletion','2 deletions','3 deletions', 'minimum 4 deletions')
  rownames(INS) <- rownames(DEL) <- id
  
  for (n in 1:length(id)) {

    
  ############# read uc
    
    
    uc <- utils::read.table(files[n])
    uc_C <- uc[which(uc$V1=='C'),]
    rownames(uc_C) <- c(1:nrow(uc_C))
    uc_C2 <- uc_C[,2]
    uc_H <- uc[which(uc$V1 =='H'),]
    rownames(uc_H) <- c(1:nrow(uc_H))
    
    
  ############# number of MSCs
    
    
    MSCs[n] <- length(uc_C2)

    
  ############# insertions
    
      
    perfectalignments <- sort(c(as.numeric(which(apply(uc_H, 1, function(x) (x[8] == paste(trimws(x[3]), 'M', sep=''))))), which(uc_H$V8 == "=")))
    PROP.PA[n] <- length(perfectalignments)/length(uc_H$V8)
      
    uc_H3 <- uc_H[-perfectalignments,]
    insertions <- gsub('[0-9]*D','', gsub('[0-9]*M','',uc_H3$V8))
    if (length(grep('I', insertions)) > 0) {
      length_insertions <- as.numeric(gsub('I','', gsub('^[I]+', '1', insertions)))
      length_insertions[is.na(length_insertions)] <- 0
      length_insertions_1bp <- length(which(length_insertions == 1))
      length_insertions_2bp <- length(which(length_insertions == 2))
      length_insertions_3bp <- length(which(length_insertions == 3))
      length_insertions_min4bp <- length(which(length_insertions >= 4))
      if(length(which(length_insertions >= 500)) > 0) {
        cat(paste0("WARNING: we found ", length(which(length_insertions >= 500)), " alignments with insertions larger than 500bp at ", id[n], "% identity. This may due to artificial minicircle dimers.\n"))
      }
        
      number_insertions <- stringr::str_count(insertions, 'I')
      number_insertions_1ins <- length(which(number_insertions == 1))
      number_insertions_2ins <- length(which(number_insertions == 2))
      number_insertions_3ins <- length(which(number_insertions == 3))
      number_insertions_min4ins <- length(which(number_insertions >= 4))
        
      INS[n,] <- c(length_insertions_1bp, length_insertions_2bp, length_insertions_3bp,length_insertions_min4bp,
          number_insertions_1ins, number_insertions_2ins, number_insertions_3ins,number_insertions_min4ins)
      } else {
        INS[n,] <- 0
      }
      

  ############# insertions
    
        
    deletions <- gsub('[0-9]*I','',gsub('[0-9]*M','',uc_H3$V8))
    if (length(grep('D', deletions)) > 0) {
      length(grep('D', deletions))
      length_deletions <- as.numeric(gsub('D','', gsub('^[D]+', '1', deletions)))
      length_deletions[is.na(length_deletions)] <- 0
      length_deletions_1bp <- length(which(length_deletions == 1))
      length_deletions_2bp <- length(which(length_deletions == 2))
      length_deletions_3bp <- length(which(length_deletions == 3))
      length_deletions_min4bp <- length(which(length_deletions >= 4))
      if(length(which(length_deletions >= 500)) > 0) {
        cat(paste0("WARNING: we found ", length(which(length_deletions >= 500)), " alignments with deletions larger than 500bp at ", id[n], "% identity. This may due to artificial minicircle dimers.\n"))
      }
        
      number_deletions <- stringr::str_count(deletions, 'D')
      number_deletions_1ins <- length(which(number_deletions == 1))
      number_deletions_2ins <- length(which(number_deletions == 2))
      number_deletions_3ins <- length(which(number_deletions == 3))
      number_deletions_min4ins <- length(which(number_deletions >= 4))
        
      DEL[n,] <- c(length_deletions_1bp, length_deletions_2bp, length_deletions_3bp,length_deletions_min4bp,
                   number_deletions_1ins, number_deletions_2ins, number_deletions_3ins,number_deletions_min4ins)
      } else {
      DEL[n,] <- 0
    }
  }
  
  names(MSCs) <- names(PROP.PA) <- id

  
  ############# plots
  
  
  theme_set(theme_minimal() + theme(axis.title = element_text(face = 'bold')))
  
  coeff <- 100/max(MSCs)
  
  gplots[[1]] <- ggplot(data.frame(MSCs), aes(x=id, y=MSCs)) +
    geom_point(alpha=0.5, size=3, color="blue") + geom_line(color="blue") +
    geom_point(data=data.frame(PROP.PA), aes(x=id, y=PROP.PA*100/coeff),alpha=0.5, size=3, color="red") + 
    geom_line(aes(x=id, y=PROP.PA*100/coeff), color="red") +
    scale_y_continuous(
      name = "Number of MSC",
      sec.axis = sec_axis(~.*coeff, name="Perfect alignment proportions (%)")
    ) + xlab("percent identity (%)") +
    theme(axis.title.y.left = element_text(color="blue"),
          axis.text.y.left = element_text(color="blue"),
          axis.title.y = element_text(color="red"),
          axis.text.y = element_text(color="red")) 
  
  melted <- reshape2::melt(INS[,1:4])
  pltINS[[1]] <- ggplot(melted, aes(x=melted$Var1, y=melted$value, color=melted$Var2)) + 
    geom_point(alpha=0.5, size=3) + geom_line() +
    xlab("percent identity (%)") + ylab("Length of insertions") + 
    theme(legend.title = element_blank())
  
  melted <- reshape2::melt(INS[,5:8])
  pltINS[[2]] <- ggplot(melted, aes(x=melted$Var1, y=melted$value, color=melted$Var2)) + 
    geom_point(alpha=0.5, size=3) + geom_line() +
    xlab("percent identity (%)") + ylab("Number of insertions") + 
    theme(legend.title = element_blank())
  
  melted <- reshape2::melt(DEL[,1:4])
  pltDEL[[1]] <- ggplot(melted, aes(x=melted$Var1, y=melted$value, color=melted$Var2)) + 
    geom_point(alpha=0.5, size=3) + geom_line() +
    xlab("percent identity") + ylab("Length of deletions") +
    theme(legend.title = element_blank())
  
  melted <- reshape2::melt(DEL[,5:8])
  pltDEL[[2]] <- ggplot(melted, aes(x=melted$Var1, y=melted$value, color=melted$Var2)) + 
    geom_point(alpha=0.5, size=3) + geom_line() +
    xlab("percent identity") + ylab("Number of deletions") +
    theme(legend.title = element_blank())
  
  gplots[[2]] <- ggpubr::ggarrange(plotlist=pltINS, ncol=2)
  gplots[[3]] <- ggpubr::ggarrange(plotlist=pltDEL, ncol=2)
  
  names(gplots) <- c("MSCs and perfect aligments", "insertions", "deletions")
  
  
  ############# return
  
  
  results <- list(MSCs, PROP.PA, INS, DEL, gplots)
  names(results) <- c("MSCs", "perfect alignments", "insertions", "deletions", "plots")
  return(results)
  
  
}



