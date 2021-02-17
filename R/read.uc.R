#' Read in uc files 
#'
#' Clustering based on a percent identity, performed with VSEARCH, generates files in uc format.
#' The read.uc function allows you to read in an uc file and store cluster and hit records in individual tables. 
#'
#' @usage read.uc(file)
#' @param file the name of the uc file e.g. all.minicircles.circ.id70.uc
#' @return
#'   \item{hits}{a table containing all hit records.}
#'   \item{clusters}{a table containing all cluster records.}
#'   \item{clustnumbers}{a vector containing the cluster numbers (0-based).}
#' @importFrom utils read.table
#' @export

read.uc <- function(file) {
  
  uc <- utils::read.table(file)
  uc_C <- uc[which(uc$V1=='C'),]
  rownames(uc_C) <- c(1:nrow(uc_C))
  uc_C2 <- uc_C[,2]
  uc_H <- uc[which(uc$V1 =='H'),]
  rownames(uc_H) <- c(1:nrow(uc_H))
  
  return(list("hits" = uc_H, "clusters" = uc_C, "clustnumbers" = uc_C2))
  
}