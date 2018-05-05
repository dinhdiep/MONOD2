library(compiler)

# Input: MHL matrix (rows: MHBs, columns: samples)
# Output; dataframe with the GSI information for each MHB (row)

gsi <- function(mhl.matrix, min.frac = 0.6) {
  mhl.matrix <- as.matrix(mhl.matrix)
  
  groups <- names(table(colnames(mhl.matrix))) # unique sample groups (i.e. brain, liver, etc.)
  index <- colnames(mhl.matrix)
  
  n.rows <- nrow(mhl.matrix)
  gsi.results <- data.frame(region = rownames(mhl.matrix), group = character(n.rows), GSI = numeric(n.rows), 
                            refMax = numeric(n.rows), stringsAsFactors = F)
  
  for (i in 1:nrow(mhl.matrix)) {
    group.means <- tapply(as.numeric(mhl.matrix[i, ]), index, function(x) mean(x, na.rm = T))
    group.max <- group.means[which.max(group.means)]
    group.assign <- names(group.max)
    
    groups.use <- groups[groups != group.assign]
    mhb.gsi <- vapply(groups.use, function(g) {
      g.mean <- mean(mhl.matrix[i, which(index == g)], na.rm = T)
      g.mhb <- 1 - (10^g.mean)/(10^group.max)
      g.mhb
    }, 1.0)
    mhb.sum <- sum(mhb.gsi, na.rm = T)/(length(groups) - 1)
    
    # Allow each MHB to be assigned to multiple tissues
    group.mean.fracs <- group.means/group.max
    if(sum(group.mean.fracs > min.frac, na.rm = T) > 1){
      group.assign <- paste(names(group.mean.fracs[group.mean.fracs > min.frac]), collapse = ",")
    }
    
    gsi.results[i, "group"] <- group.assign
    gsi.results[i, "GSI"] <- mhb.sum
    gsi.results[i, "refMax"] <- group.max
  }
  
  return(gsi.results)
}
gsi <- cmpfun(gsi) # Compile gsi function for faster runtime


pdf("results/plasma_prediction/plots.pdf")
wgbs.meta.file = "data/ng.3805/WGBS.getHaplo.sampleInfo.txt"; 
wgbs.meta <- read.table(wgbs.meta.file, sep ='\t', header=T, stringsAsFactors = F)

wgbs.normal_samples.meta <- wgbs.meta[(wgbs.meta$Disease.Type == "normal"),]
wgbs.normal_samples.meta$Tissue.Type <- gsub("small intestine", "GI", wgbs.normal_samples.meta$Tissue.Type)
wgbs.normal_samples.meta$Tissue.Type <- gsub("colon", "GI", wgbs.normal_samples.meta$Tissue.Type)
wgbs.normal_samples.meta$Tissue.Type <- gsub("esophagus", "GI", wgbs.normal_samples.meta$Tissue.Type)
wgbs.normal_samples.meta$Tissue.Type <- gsub("stomach", "GI", wgbs.normal_samples.meta$Tissue.Type)

normal_blood_vessel.meta <- wgbs.normal_samples.meta[(wgbs.normal_samples.meta$Tissue.Type == "blood")|(wgbs.normal_samples.meta$Tissue.Type == "vessel"),]
normal_blood.meta <- wgbs.normal_samples.meta[(wgbs.normal_samples.meta$Tissue.Type == "blood"),]

exclude.samples <- c("SAMN00847541")

keep.tissues <- c("GI", "neural", "heart", "liver", "pancreas", "kidney", "fat", "lung", "pancreas")
keep.tissues.meta <- wgbs.normal_samples.meta[((wgbs.normal_samples.meta$Tissue.Type %in% keep.tissues) & !(wgbs.normal_samples.meta$Sample.ID %in% exclude.samples)),]
cancer_samples.meta <- wgbs.meta[(((wgbs.meta$Disease.Type == "cancer")|(wgbs.meta$Disease.Type == "cancer metastasis"))&(!(wgbs.meta$Sample.ID %in% exclude.samples))),]

repetitive_MHB_file <- "data/ng.3805/overlap_rmsk.mhbsv1.txt"
repetitive_MHB_list <- read.table(repetitive_MHB_file, sep='\t', header=F, stringsAsFactors = F)

wgbs.umhl.file = "data/ng.3805/WGBS_180118.gethaplo.umhl.mhbs1.0.useSampleID.txt.gz"; 
wgbs.umhl.dat <- read.table(gzfile(wgbs.umhl.file), sep='\t', header=T, row.names=1)
dim(wgbs.umhl.dat)
wgbs.mhl.file = "data/ng.3805/WGBS_180118.gethaplo.mhl.mhbs1.0.useSampleID.txt.gz"; 
wgbs.mhl.dat <- read.table(gzfile(wgbs.mhl.file), sep='\t', header=T, row.names=1)
dim(wgbs.mhl.dat)

wgbs.umhl.dat <- wgbs.umhl.dat[!(rownames(wgbs.umhl.dat) %in% repetitive_MHB_list$V1),]
dim(wgbs.umhl.dat)
rownames(wgbs.umhl.dat) <- paste0("uMHL:",rownames(wgbs.umhl.dat))

wgbs.mhl.dat <- wgbs.mhl.dat[!(rownames(wgbs.mhl.dat) %in% repetitive_MHB_list$V1),]
dim(wgbs.mhl.dat)
rownames(wgbs.mhl.dat) <- paste0("MHL:",rownames(wgbs.mhl.dat))

wgbs.combined.dat <- rbind(wgbs.mhl.dat,wgbs.umhl.dat)
all_included_samples <- c(keep.tissues.meta$Sample.ID, normal_blood_vessel.meta$Sample.ID, cancer_samples.meta$Sample.ID)
print(all_included_samples)
wgbs.combined.dat <- wgbs.combined.dat[,all_included_samples]
dim(wgbs.combined.dat)

wgbs.combined.dat.marker.nsamples <- apply(wgbs.combined.dat, 1, function(x) sum(!is.na(x)))
wgbs.combined.dat.sample.nmarkers <- apply(wgbs.combined.dat, 2, function(x) sum(!is.na(x)))

hist(wgbs.combined.dat.marker.nsamples)
hist(wgbs.combined.dat.sample.nmarkers)

wgbs.filtered.combined.dat = wgbs.combined.dat[(wgbs.combined.dat.marker.nsamples>58),]
dim(wgbs.filtered.combined.dat)

keep.tissues.combined.dat <- wgbs.filtered.combined.dat[,keep.tissues.meta$Sample.ID]
dim(keep.tissues.combined.dat)
keep.tissues.markers_ninty_pct <- apply(keep.tissues.combined.dat , 1, function(x) quantile(x, 0.9, na.rm = TRUE))
hist(keep.tissues.markers_ninty_pct)

normal_blood.combined.dat <- wgbs.filtered.combined.dat[,normal_blood.meta$Sample.ID]
dim(normal_blood.combined.dat)
normal_blood.markers_ninety_pct <- apply(normal_blood.combined.dat, 1, function(x) quantile(x, 0.9, na.rm = TRUE))
hist(normal_blood.markers_ninety_pct)

cancer_samples.combined.dat <- wgbs.filtered.combined.dat[,cancer_samples.meta$Sample.ID]
dim(cancer_samples.combined.dat)
cancer_samples.markers_fifty_pct <- apply(cancer_samples.combined.dat , 1, function(x) quantile(x, 0.5, na.rm = TRUE))
hist(cancer_samples.markers_fifty_pct)

keep.tissues.low_in_blood_vessel.some_in_other.combined.dat <-wgbs.filtered.combined.dat[((normal_blood.markers_ninety_pct<0.2)&((cancer_samples.markers_fifty_pct > 0.3)|(keep.tissues.markers_ninty_pct > 0.3))),]
dim(keep.tissues.low_in_blood_vessel.some_in_other.combined.dat)

colnames(keep.tissues.low_in_blood_vessel.some_in_other.combined.dat) <- c(keep.tissues.meta$Tissue.Type, rep("blood", length(normal_blood_vessel.meta$Sample.ID)), rep("cancer", length(cancer_samples.meta$Sample.ID)))
print(colnames(keep.tissues.low_in_blood_vessel.some_in_other.combined.dat))
print(colnames(wgbs.filtered.combined.dat))
keep.tissues.low_in_blood.feature.gsi <- gsi(keep.tissues.low_in_blood_vessel.some_in_other.combined.dat, min.frac = 0.8) # Calculate GSI for each feature
keep.tissues.low_in_blood.feature.gsi <- keep.tissues.low_in_blood.feature.gsi[order(keep.tissues.low_in_blood.feature.gsi$GSI, decreasing = T),]

write.table(keep.tissues.low_in_blood.feature.gsi, "results/plasma_prediction/WGBS_MHBv1_MHL_low_in_blood_9-tissue_with_cancer_GSI.list.txt", sep = '\t', col.names = T, row.names = F)
dev.off()

