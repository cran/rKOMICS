#' Cluster Analyses
#' 
#' The function msc.uc reads the output of the clustering analyses (UC file) for each specified minimum percent identity (MPI) into a single list, which will be analyzed automatically to calculate and visualize, per MPI, the number of minicircle sequence classes (MSCs), the proportion of perfect alignments (i.e. alignments without any insertion/deletion, but allowing point mutations) and the number of alignment gaps. Gaps are defined by i) the number of insertions/deletions and ii) the length in base pairs of each individual insertion/deletion. It also issues a warning when large gaps (>500 bp) are found, which points the user to anomalous alignments due to e.g. artificial dimers introduced by the assembly process. This allows the user to make an informed decision about the MPI (or MPI's) that best captures minicircle sequence richness within a (group of) sample(s) while minimizing the number and length of alignment gaps.
#' 
#' @usage msc.uc(files)
#' @param files a character vector that includes the file names of UC files (produced by USEARCH or VSEARCH), such as all.minicircles.circ.id70.uc, all.minicircles.circ.id80.uc, and so on. Please ensure that your file names end with 'idxx.uc' for this function to work properly.
#' @return
#'   \item{MSCs}{a numerical vector containing the number of MSC per MPI.}
#'   \item{perfect aligments}{a numerical vector containing the proportions of perfect alignments per MPI.}
#'   \item{insertions}{a list showing the insertion lengths per MPI. Each element in the list corresponds to a specific MPI, and it provides the lengths of identified insertions.}
#'   \item{deletions}{a list showing the deletions lengths per MPI. Each element in the list corresponds to a specific MPI, and it provides the lengths of identified deletions.}
#'   \item{insertions summary}{a table showing the length and the number of insertions across different MPIs.}
#'   \item{deletions summary}{a table showing the length and the number of deletions across different MPIs.}
#'   \item{plots}{various plots showing previous results.}
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


  #############   0   Tests   #############

  
  if (sum(file.exists(files))!=length(files)) 
    stop(paste0("One or more files don't exist in ", getwd()))
  if (length(files)==0) 
    stop("files argument is empty")
  if (length(grep("id", files)) != length(files)) 
    stop("Please ensure that your file names end with 'idxx.uc' for this function to work properly.")
  
  
  #############   0   Define global functions or variables   #############
  
  
  Var1 <- value <- Var2 <- NULL
  
  
  #############   1   Analysis   #############
  
  
  ids <- as.numeric(gsub("\\D", "", tools::file_path_sans_ext(files)))
  
  MSCs <- PROP.PA <- vector("numeric", length(ids))
  gplots <- pltINS <- pltDEL <- insertionsA <- deletionsA <- list()
  INS <- DEL <- matrix(ncol = 8, nrow = length(ids))
  colnames(INS) <- c('1bp','2bp','3bp','minimum 4bp',
                     '1 insertion','2 insertions','3 insertions', 'minimum 4 insertions')
  colnames(DEL) <- c('1bp','2bp','3bp','minimum 4bp',
                     '1 deletion','2 deletions','3 deletions', 'minimum 4 deletions')
  rownames(INS) <- rownames(DEL) <- ids
  
  for (n in seq_along(ids)) {

    
  #############   2   Read uc   #############
    
    
    uc <- utils::read.table(files[n])
    uc_C <- uc[uc$V1 == 'C', , drop = FALSE]
    rownames(uc_C) <- seq_len(nrow(uc_C))
    uc_C2 <- uc_C[, 2]
    uc_H <- uc[uc$V1 == 'H', , drop = FALSE]
    rownames(uc_H) <- seq_len(nrow(uc_H))
    
    
  #############   3   Number of MSCs   #############
    
    
    MSCs[n] <- length(uc_C2)

    
  #############   4   Perfect alignments    #############
    
    
    pa <- uc_H$V8 == paste(trimws(uc_H$V3), 'M', sep = '')
    perfectalignments <- sort(c(which(pa), which(uc_H$V8 == "=")))
    PROP.PA[n] <- length(perfectalignments) / length(uc_H$V8)
  

  #############   5   Insertion   #############
    
    
    uc_H3 <- uc_H[-perfectalignments, , drop = FALSE]
    insertions <- gsub('[0-9]*D','', gsub('[0-9]*M','',uc_H3$V8))
    
    if (any(grepl('I', insertions))) {
      insertions_split <- unlist(stringr::str_extract_all(insertions, '\\d+[I]|\\D[I]|[I]'))
      insertionsA[[n]] <- as.numeric(gsub('I', '', insertions_split[grep('\\d',insertions_split)]))
      insertions_short <- stringr::str_count(insertions_split[!grepl("\\d", insertions_split)], "I")
      insertions_1bp <- length(insertions_short)
      insertions_2bp <- sum(insertionsA[[n]] == 2)
      insertions_3bp <- sum(insertionsA[[n]] == 3)
      insertions_min4bp <- sum(insertionsA[[n]] >= 4)
      
      if (sum(insertionsA[[n]] >= 500) > 0) {
        cat(paste0("WARNING: we found ", sum(insertionsA[[n]] >= 500), 
                   " alignments with insertions larger than 500bp at ", 
                   ids[n], 
                   "% identity. This may due to artificial minicircle dimers.\n"))
      }
      
      insertions_short <- c(insertions_short, rep(1,length(insertions_short[insertions_short==2])))
      insertions_short[insertions_short==2] <- 1
      insertionsA[[n]] <- sort(c(insertionsA[[n]], insertions_short))
      
      number_insertions <- stringr::str_count(insertions, 'I')
      number_insertions_1ins <- sum(number_insertions == 1)
      number_insertions_2ins <- sum(number_insertions == 2)
      number_insertions_3ins <- sum(number_insertions == 3)
      number_insertions_min4ins <- sum(number_insertions >= 4)
      
      INS[n, ] <- c(
        insertions_1bp,
        insertions_2bp,
        insertions_3bp,
        insertions_min4bp,
        number_insertions_1ins,
        number_insertions_2ins,
        number_insertions_3ins,
        number_insertions_min4ins
      )
      
    } else {
      INS[n,] <- 0
      insertionsA[[n]] <- 0
    }

    
  #############   6   Deletions   #############
    
        
    deletions <- gsub('[0-9]*I','',gsub('[0-9]*M','',uc_H3$V8))
    
    if (any(grepl('D', deletions))) {
      deletions_split <- unlist(stringr::str_extract_all(deletions, '\\d+[D]|\\D[D]|[D]'))
      deletionsA[[n]] <- as.numeric(gsub('D', '', deletions_split[grep('\\d',deletions_split)]))
      deletions_short <- stringr::str_count(deletions_split[!grepl("\\d", deletions_split)], "D")
      deletions_1bp <- length(deletions_short)
      deletions_2bp <- sum(deletionsA[[n]] == 2)
      deletions_3bp <- sum(deletionsA[[n]] == 3)
      deletions_min4bp <-  sum(deletionsA[[n]] >= 4)
      
      if (sum(deletionsA[[n]] >= 500) > 0) {
        cat(paste0("WARNING: we found ", 
                   sum(deletionsA[[n]] >= 500), 
                   " alignments with deletions larger than 500bp at ", 
                   ids[n], 
                   "% identity. This may be due to artificial minicircle dimers.\n"))
      }
      
      deletions_short <- c(deletions_short, rep(1,length(deletions_short[deletions_short==2])))
      deletions_short[deletions_short==2] <- 1
      deletionsA[[n]] <- sort(c(deletionsA[[n]], deletions_short))
      
      number_deletions <- stringr::str_count(deletions, 'D')
      number_deletions_1ins <- sum(number_deletions == 1)
      number_deletions_2ins <- sum(number_deletions == 2)
      number_deletions_3ins <- sum(number_deletions == 3)
      number_deletions_min4ins <- sum(number_deletions >= 4)
      
      DEL[n, ] <- c(
        deletions_1bp,
        deletions_2bp,
        deletions_3bp,
        deletions_min4bp,
        number_deletions_1ins,
        number_deletions_2ins,
        number_deletions_3ins,
        number_deletions_min4ins
      )
      
    } else {
      DEL[n,] <- 0
      deletionsA[[n]] <- 0
    }
    
  }
  
  names(MSCs) <- names(PROP.PA) <- names(insertionsA) <- names(deletionsA) <- ids
    
    
    #############   7   Plots   #############
  
  
  theme_set(theme_minimal() + 
              theme(axis.title = element_text(face = 'bold'),
                    legend.title = element_blank()))
  
  coeff <- 100/max(MSCs)
  
  gplots[[1]] <- ggplot(data.frame(MSCs), aes(x = ids, y = MSCs)) +
    geom_point(alpha = 0.8, size = 3, color = "blue") +
    geom_line(color = "blue") +
    geom_point(
      data = data.frame(PROP.PA),
      aes(x = ids, y = PROP.PA * 100 / coeff),
      alpha = 0.5,
      size = 3,
      color = "red"
    ) +
    geom_line(aes(x = ids, y = PROP.PA * 100 / coeff), color = "red") +
    scale_y_continuous(
      name = "Number of MSC",
      sec.axis = sec_axis(~. * coeff, name = "Perfect alignment proportions (%)")
    ) +
    xlab("percent identity (%)") +
    theme(
      axis.title.y.left = element_text(color = "blue"),
      axis.text.y.left = element_text(color = "blue"),
      axis.title.y = element_text(color = "red"),
      axis.text.y = element_text(color = "red")
    )
  
  
  pltINS[[1]] <- ggplot(na.omit(reshape2::melt(INS[, 1:4])), aes(x = Var1, y = value, color = Var2)) +
    geom_point(alpha = 0.8, size = 3) +
    geom_line() +
    xlab("percent identity (%)") +
    ylab("Total number of insertions") +
    theme(legend.position = 'bottom')
  
  pltINS[[2]] <- ggplot(na.omit(reshape2::melt(INS[, 5:8])), aes(x = Var1, y = value, color = Var2)) +
    geom_point(alpha = 0.8, size = 3) +
    geom_line() +
    xlab("percent identity (%)") +
    ylab("Number of alignments") +
    theme(legend.position = 'bottom')
  
  pltDEL[[1]] <- ggplot(na.omit(reshape2::melt(DEL[, 1:4])), aes(x = Var1, y = value, color = Var2)) +
    geom_point(alpha = 0.8, size = 3) +
    geom_line() +
    xlab("percent identity (%)") +
    ylab("Total number of deletions") +
    theme(legend.position = 'bottom')
  
  pltDEL[[2]] <- ggplot(na.omit(reshape2::melt(DEL[, 5:8])), aes(x = Var1, y = value, color = Var2)) +
    geom_point(alpha = 0.8, size = 3) +
    geom_line() +
    xlab("percent identity (%)") +
    ylab("Number of alignments") +
    theme(legend.position = 'bottom')
  
  gplots[[2]] <- ggpubr::ggarrange(plotlist = pltINS, ncol = 2)
  gplots[[3]] <- ggpubr::ggarrange(plotlist = pltDEL, ncol = 2)
  
  names(gplots) <- c("MSCs and perfect aligments", "insertions", "deletions")
  
  
  #############   8   Return   #############
  
  
  results <- list(insertionsA, deletionsA, 
                  MSCs, PROP.PA, INS, DEL, gplots)
  names(results) <- c("insertions", "deletions", 
                      "MSCs", "perfect alignments", 
                      "insertions summary", "deletions summary", "plots")
  return(results)
  
  
}



