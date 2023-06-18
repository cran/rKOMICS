#' Read in uc files 
#'
#' The read.uc function is used to read the output of clustering analyses from a UC file. The function stores the information in two tables: one for hit records (H) and one for cluster records (C).
#'
#' @usage read.uc(file)
#' @param file the name of the UC file that contains the clustering analysis results. The file should be specified with its full name, including the extension (e.g., all.minicircles.circ.id70.uc).
#' @return
#'   \item{hits}{a table containing all the hit records from the UC file. Each row of the table represents a hit record, providing information about the alignment between a query sequence and a target sequence.}
#'   \item{clusters}{a table containing all the cluster records from the UC file. Each row of the table represents a cluster record, providing information about the clustering of sequences into clusters.}
#'   \item{clustnumbers}{a vector containing the cluster numbers (0-based). Each element of the vector represents a cluster identified in the clustering analysis.}
#' @importFrom utils read.table
#' @export

read.uc <- function(file) {
  
  uc <- utils::read.table(file)
  uc_C <- uc[which(uc$V1=='C'),]
  rownames(uc_C) <- c(1:nrow(uc_C))
  uc_C2 <- uc_C[,2]
  uc_H <- uc[which(uc$V1 =='H'),]
  if (nrow(uc_H) != 0) {
    rownames(uc_H) <- c(1:nrow(uc_H))
  }
  
  return(list("hits" = uc_H, "clusters" = uc_C, "clustnumbers" = uc_C2))
  
}
