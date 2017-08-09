library(MASS)
library(pheatmap)
library(RColorBrewer)
library(randomForest)
library(ggplot2)
library(impute)
library(earth)
library(ModelMetrics)
library(abind)


# Note: the below analysis is performed with OUTLIERS <- c('.NONE.') and
# outliers observed on PC1/PC2 get added to the outlier list; and then it
# is re-run.

OUTLIERS <- c('CRC.T.4.2014')

args <- commandArgs(trailingOnly = TRUE)
myrrbs.file <- args[1]
mywgbs.file <- args[2]
myrrbs.meta.file <- args[3]
mywgbs.meta.file <- args[4]
name.rdata <- args[5]

######
# LOADING DATA
######

gen.data <- function(rrbs.meta.file = "RRBS.getHaplo.sampleInfo.txt", 
                     wgbs.meta.file = "WGBS.getHaplo.sampleInfo.txt",
			               rrbs.file = "RRBS_170609.gethaplo.mhl.mhbs1.0.useSampleID.txt", 
			               wgbs.file = "WGBS.getHaplo.mhl.mhbs1.0.rmdup_consistent.useSampleID.txt",
			               MIN.N.OBSERVED.SAMPLES = 50000,  MIN.N.OBSERVED.SITES = 60,
			               MIN.SAMPLE.NMARKERS = 80000, MIN.MARKER.NSAMPLES = 70) {
  
  rrbs.meta <- read.table(rrbs.meta.file, header = T, sep='\t', stringsAsFactors = F)
  wgbs.meta <- read.table(wgbs.meta.file, sep ='\t', header=T, stringsAsFactors = F)
  
  rrbs.dat <- read.table(gzfile(rrbs.file), sep='\t', header=T, row.names=1)
  wgbs.dat <- read.table(gzfile(wgbs.file), sep='\t', header=T, row.names=1)
  
  rrbs.meta$Cancer.Origin <- rrbs.meta$Tissue ## rename this
  
  rownames(wgbs.meta) <- make.names(wgbs.meta$Sample.ID)
  
  intersect.names <- sort(intersect(colnames(wgbs.dat), rownames(wgbs.meta)))
  
  wgbs.dat <- wgbs.dat[,intersect.names]
  wgbs.meta <- wgbs.meta[intersect.names,]
  
  rrbs.meta <- subset(rrbs.meta, ! rrbs.meta$Sample %in% OUTLIERS)
  
  rrbs.matched <- rrbs.dat[,make.names(rrbs.meta$Sample)]
  
  rrbs.observed.marker <- apply(rrbs.matched, 1, function(x) { sum(! is.na(x))})
  rrbs.meta$observed.sites <- apply(rrbs.matched, 2, function(x) { sum(! is.na(x))})
  pdf('20170609.observed.sites.per.sample.pdf')
  gp <- ggplot(rrbs.meta, aes(x=Disease, y=observed.sites, fill=Tissue, color=Type)) + geom_boxplot()
  print(gp)
  dev.off()
  
  MIN.N.OBSERVED.SAMPLES <- 50000
  rrbs.observed.marker <- data.frame('observed'=rrbs.observed.marker, 
                                     'chrom'=sapply(strsplit(names(rrbs.observed.marker), ':'),  function(x) x[1]))
  
  pdf('20170609.site.coverage.pdf')
  gp <- ggplot(rrbs.observed.marker, aes(x=chrom, y=observed, fill=chrom)) + geom_boxplot()
  print(gp)
  dev.off()
  
  MIN.N.OBSERVED.SITES <- 60
  
  min.obs.samp <- min(rrbs.meta$observed.sites)
  min.obs.site <- min(rrbs.observed.marker$observed)
  
  rrbs.observed.marker <- rrbs.observed.marker$observed
  rrbs.meta.sub <- rrbs.meta
  rrbs.matched.sub <- rrbs.matched
  
  
  #############
  # RRBS CLEANING
  #############
  
  # iteratively exclude samples and sites until all samples and sites cross the threshold
  while ( min.obs.samp < MIN.N.OBSERVED.SAMPLES && min.obs.site < MIN.N.OBSERVED.SITES) {
    rrbs.meta.sub <- rrbs.meta.sub[rrbs.meta.sub$observed.sites >= MIN.N.OBSERVED.SAMPLES,]
    rrbs.matched.sub <- rrbs.matched.sub[rrbs.observed.marker>= MIN.N.OBSERVED.SITES, make.names(rrbs.meta.sub$Sample)]
    rrbs.observed.marker <- apply(rrbs.matched.sub, 1, function(x) { sum(! is.na(x))})
    rrbs.meta.sub$observed.sites <- apply(rrbs.matched.sub, 2, function(x) { sum(! is.na(x))})
    
    min.obs.samp <- min(rrbs.meta.sub$observed.sites)
    min.obs.site <- min(rrbs.observed.marker)
    
  }
  
  # now re-plot
  
  pdf('20170609.postfilter.observed.sites.per.sample.pdf')
  gp <- ggplot(rrbs.meta.sub, aes(x=Disease, y=observed.sites, fill=Tissue, color=Type)) + geom_boxplot()
  print(gp)
  dev.off()
  
  pdf('20170609.postfilter.site.coverage.pdf')
  rrbs.observed.marker <- data.frame('observed'=rrbs.observed.marker, 
                                     'chrom'=sapply(strsplit(names(rrbs.observed.marker), ':'),  function(x) x[1]))
  gp <- ggplot(rrbs.observed.marker, aes(x=chrom, y=observed, fill=chrom)) + geom_boxplot()
  print(gp)
  dev.off()
  
  
  #################################
  ### WGBS CLEANING
  ################################
  
  wgbs.marker.nsamples <- apply(wgbs.dat, 1, function(x) { sum(! is.na(x))})
  wgbs.sample.nmarkers <- apply(wgbs.dat, 2, function(x) { sum(! is.na(x))})
  
  hist(wgbs.sample.nmarkers, 30)
  hist(wgbs.marker.nsamples, 300)
  
  ## iteratively exclude markers and samples until the above thresholds are met
  wgbs.dat.sub <- wgbs.dat
  while(any(wgbs.marker.nsamples < MIN.MARKER.NSAMPLES) || any(wgbs.sample.nmarkers < MIN.SAMPLE.NMARKERS) ) {
    wgbs.dat.sub <- wgbs.dat.sub[wgbs.marker.nsamples >= MIN.MARKER.NSAMPLES, wgbs.sample.nmarkers >= MIN.SAMPLE.NMARKERS]
    wgbs.marker.nsamples <- apply(wgbs.dat.sub, 1, function(x) sum(!is.na(x)))
    wgbs.sample.nmarkers <- apply(wgbs.dat.sub, 2, function(x) sum(!is.na(x)))
  }
  
  wgbs.meta.sub <- wgbs.meta[colnames(wgbs.dat.sub),]
  print(dim(wgbs.meta.sub))
  
  table(wgbs.meta.sub$Tissue.Type)
  # place `bad` tissues (n <= 3) into `other` category
  low.n.tissues <- names(table(wgbs.meta.sub$Tissue.Type))[table(wgbs.meta.sub$Tissue.Type) < 4]
  wgbs.meta.sub$Clean.Tissues <- sapply(wgbs.meta.sub$Tissue.Type, function(x) { if ( x %in% low.n.tissues ) 'Other' else x})
  table(wgbs.meta.sub$Clean.Tissues)
  
  ## subset to RRBS markers
  wgbs.dat.sub <- wgbs.dat.sub[rownames(rrbs.matched.sub)[rownames(rrbs.matched.sub) %in% rownames(wgbs.dat.sub)],]
  
  
  #####################
  ### WGBS IMPUTATION (K-NN)
  #####################
  
  scale.sdf <- function(X) {
    U <- scale(X)
    isnan <- which(apply(U, 2, function(x) { any(is.na(x))}))
    if ( length(isnan) > 0 ) {
      U <- U[,-isnan]
    }
    
    U
  }
  
  wgbs.imp <- impute.knn(as.matrix(wgbs.dat.sub))
  wgbs.imp <- wgbs.imp$data
  print(dim(wgbs.imp))
  
  ######################
  ### RRBS IMPUTATION (K-NN)
  ######################
  set.seed(1618033) # golden ratio
  rrbs.matched.imp <- impute.knn(as.matrix(rrbs.matched.sub), k=10, rowmax=0.8)
  rrbs.matched.imp <- rrbs.matched.imp$data
  colnames(rrbs.matched.imp) <- colnames(rrbs.matched.sub)
  rownames(rrbs.matched.imp) <- rownames(rrbs.matched.sub)
 
  print(dim(rrbs.matched.imp)) 
  
  ######################
  ### DEFINE INTERSECTIONAL DATA
  ######################
  
  rrbs.geneset <- rownames(rrbs.matched.imp[rrbs.observed.marker$observed > 0.7*nrow(rrbs.meta.sub), rrbs.meta.sub$Type == 'Plasma'])
  wgbs.geneset <- intersect(rrbs.geneset, rownames(wgbs.imp))
  wgbs.imp.pred <- wgbs.imp[wgbs.geneset,]
  rrbs.imp.pred <- rrbs.matched.imp[rrbs.geneset, ]

  list(rrbs.meta=rrbs.meta.sub, rrbs.dat=rrbs.imp.pred, rrbs.matched.imp=rrbs.matched.imp, wgbs.meta=wgbs.meta.sub, wgbs.dat=wgbs.imp.pred)
}

#### Run to generate imputed matrices ####

orig.data <- gen.data(rrbs.file = myrrbs.file,
                      wgbs.file = mywgbs.file,
                      rrbs.meta.file = myrrbs.meta.file, 
                      wgbs.meta.file = mywgbs.meta.file)

save.image(name.rdata)
