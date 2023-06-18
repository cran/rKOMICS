### R code from vignette source 'example.Rnw'

###################################################
### code chunk number 1: example.Rnw:45-48
###################################################
library(rKOMICS)
data(exData, package = "rKOMICS")
table(exData$species)


###################################################
### code chunk number 2: example.Rnw:65-68
###################################################
library(ggplot2)
library(ggpubr)
library(viridis)


###################################################
### code chunk number 3: example.Rnw:89-91
###################################################
map <- msc.quality(mapstats = system.file("extdata", exData$mapstats, package = "rKOMICS"),exData$species)
lapply(map$proportions, mean)$MR_HQ 


###################################################
### code chunk number 4: example.Rnw:98-99
###################################################
map$plots$MR_HQ + labs(caption = paste0('Proportion of mapped reads with high quality, ', Sys.Date()))


###################################################
### code chunk number 5: example.Rnw:123-129
###################################################
depth <- msc.depth(depthstats = system.file("extdata", 
                                            exData$depthstats, 
                                            package = "rKOMICS"), 
                   groups = exData$species,
                   HCN = exData$medGWD/2) # haploid copy number
depth$CN


###################################################
### code chunk number 6: example.Rnw:150-152
###################################################
bf <- msc.length(file = system.file("extdata", "all.minicircles.fasta", package = "rKOMICS"), samples = exData$samples,groups = exData$subspecies)
bf$plot 


###################################################
### code chunk number 7: example.Rnw:162-163
###################################################
c(length(bf$length),length(which(bf$length < 800)),length(which(bf$length > 1400)))


###################################################
### code chunk number 8: example.Rnw:185-187
###################################################
pre <- preprocess(files = system.file("extdata", exData$fastafiles,package = "rKOMICS"),groups = exData$species,circ = TRUE, min = 500, max = 1200, writeDNA = FALSE)
pre$summary


###################################################
### code chunk number 9: example.Rnw:195-196
###################################################
pre$plot + labs(caption = paste0('N of MC sequences before and after filtering, ', Sys.Date()))


###################################################
### code chunk number 10: example.Rnw:214-215
###################################################
ucs <- msc.uc(files = system.file("extdata", exData$ucs, package = "rKOMICS"))


###################################################
### code chunk number 11: example.Rnw:218-219
###################################################
c(ucs$MSCs["100"],ucs$MSCs["97"],ucs$MSCs["95"])


###################################################
### code chunk number 12: example.Rnw:226-227
###################################################
ucs$plots$`MSCs and perfect aligments`


###################################################
### code chunk number 13: example.Rnw:238-239
###################################################
ucs$plots$insertions


###################################################
### code chunk number 14: example.Rnw:263-265
###################################################
data(matrices, package = "rKOMICS")
colSums(matrices[["id97"]]) # --> number of MSC per sample


###################################################
### code chunk number 15: example.Rnw:276-279
###################################################
msc.heatmap(clustmatrix = matrices[["id97"]], 
            groups = exData$species,
            samples = exData$samples)


###################################################
### code chunk number 16: example.Rnw:298-300
###################################################
richness <- msc.richness(matrices, samples = exData$samples, groups = exData$species)
apply(richness$table[which(richness$table$group=="L. peruviana"),-(1:2)], 2, mean)


###################################################
### code chunk number 17: example.Rnw:307-308
###################################################
apply(richness$table[which(richness$table$group=="L. braziliensis"),-(1:2)], 2, mean)


###################################################
### code chunk number 18: example.Rnw:315-316
###################################################
apply(richness$table[which(richness$table$group=="hybrid"),-(1:2)], 2, mean) 


###################################################
### code chunk number 19: example.Rnw:323-324
###################################################
richness$plot


###################################################
### code chunk number 20: example.Rnw:344-347
###################################################
sim <- msc.similarity(matrices, samples = exData$samples, 
                      groups = exData$species)
sim$relfreq.plot + scale_fill_viridis(discrete = TRUE)


###################################################
### code chunk number 21: example.Rnw:356-358
###################################################
c(sim$relfreq$id97["2"]*100,
  sim$relfreq$id97["3"]*100)


###################################################
### code chunk number 22: example.Rnw:376-379
###################################################
res.pca <- lapply(matrices, function(x) msc.pca(x, samples = exData$samples, 
                  groups = exData$species, n=30, labels=FALSE, title=NULL))
res.pca$id97$plot


