#' Build cluster matrix
#'
#' Clustering based on a percent identity, performed with the VSEARCH tool, will generate files in uc format.
#' The msc.matrix function will transform every input file into a cluster matrix.
#' The columns of the matrix correspond to the samples and the rows of the matrix correspond to the minicircle sequence cluster (MSC).
#' The absence of a MSC in a sample is indicated with the value of zero while the presence of a MSC in a sample will be indicated with a value >= 1.
#'
#' @usage msc.matrix(files, samples, groups)
#' @param files a character vector containing the uc file names (output of the VSEARCH tool) e.g. all.minicircles.circ.id70.uc, all.minicircles.circ.id80.uc...
#' @param samples a character vector containing the sample names. The vector should be order alphabetically.
#' @param groups a vector, of equal length as samples, specifying to which group (e.g. species) the samples belong to.
#' @return a list containing one cluster matrix per percent identity. The value 0 in the cluster matrix means the MSC doesn't occur in the sample. A value higher than 0 means the MSC is found at least once in the sample.
#' @examples
#' data(exData)
#'
#' ### run function
#' \donttest{
#' matrices <- msc.matrix(files = system.file("extdata", exData$ucs, package="rKOMICS"), 
#'                       samples = sort(exData$samples), 
#'                       groups = exData$species[order(exData$samples)])
#' }
#' 
#' ### or: 
#' data(matrices)
#' 
#' ### show matrix with id 95%
#' matrices[["id95"]]
#' rowSums(matrices[["id95"]]) # --> frequency of MSC across all samples
#' colSums(matrices[["id95"]]) # --> number of MSC per sample
#'
#' @export

msc.matrix <- function(files, samples, groups) {


  ############# tests

  
  if (sum(file.exists(files))!=length(files)) stop("ERROR: One or more files doesn't exist")
  if (!identical(samples, sort(samples))) stop("ERROR: Make sure to order you sample vector. Adapt your groups factor if necessary.")
  if (length(samples) != length(groups)) stop("ERROR: The sample and groups vector are not of equal length")
  if (!is.factor(groups)) stop("ERROR: The groups vector should be a factor")


  ############# create cluster matrix


  id <- as.numeric(regmatches( gsub(".*/", "", files), gregexpr("[[:digit:]]+",  gsub(".*/", "", files))))
  
  clust_mat_ids <- list()

  for (n in 1:length(id)) {
    
    cat(paste0("Calculating cluster matrix with id", id[n], "...\n"))
    
    uc_H <- read.uc(files[n])$hits
    uc_C <- read.uc(files[n])$clusters
    uc_C2 <- read.uc(files[n])$clustnumbers

    clusters <- vector(mode = 'list', length = length(uc_C))
    for (i in 1:length(uc_C2)) {
      if (uc_C[uc_C$V2==uc_C2[i], 3] == 1) {
        clusters[[i]] <- as.character(uc_C[uc_C$V2==uc_C2[i], 9])
      } else {
        clusters[[i]] <- unique(c(as.character(uc_H[uc_H$V2==uc_C2[i], 9]), as.character(uc_H[uc_H$V2==uc_C2[i],10])))
      }
    }

    clust_mat <- matrix(nrow=length(clusters), ncol = length(samples))
    colnames(clust_mat) <- samples
    rownames(clust_mat) <- paste('C', uc_C[,2], sep = '')
    for (t in 1:length(samples)) {
      temp <- lapply(lapply(clusters, function(x) table(gsub('_contig.*','',x))),
                     function(x) x[which(names(x) == samples[t])])

      for (ct in 1:length(temp)) {
        if (length(temp[[ct]]) > 0) {
          clust_mat[ct,t] <- temp[[ct]]
        } else if (length(temp[[ct]]) == 0) {
          clust_mat[ct,t] <- 0 #temp[[ct]]
        }
      }
    }

    clust_mat_ids[[n]] <- clust_mat[,order(colnames(clust_mat))]
  }

  names(clust_mat_ids) <- paste0("id", id)
  return(clust_mat_ids)

}
