library(caret)
library(ggplot2)
library(pROC)

result_dir <- "results/plasma_prediction"
data_dir <- "data/ng.3805"

pdf(paste0(result_dir, "/plasma_prediction_plots.pdf"), width = 3, height = 3)

seed.table <- read.table(paste0(data_dir, "/seed.table.txt"), header = F)$V1[1:1]

normal.vs.cancer.train <- read.table(paste0(data_dir, "/Normal_vs_cancer_train.cancer.samples.txt"), header = T)$x
normal.vs.cancer.test <- read.table(paste0(data_dir, "/Normal_vs_cancer_test.cancer.samples.txt"), header = T)$x
colon.vs.lung.train <- read.table(paste0(data_dir, "/Colon_vs_lung_train.types.samples.txt"), header = T)$x
colon.vs.lung.test <- read.table(paste0(data_dir, "/Colon_vs_lung_test.types.samples.txt"), header = T)$x

args <- commandArgs(trailingOnly = TRUE)

OUTLIERS <- c( 'NC.P.1.2014', 
          'NC.P.20.2014', 
          'NC.P.21.2014', 
          'NC.P.24.2014', 
          'NC.P.3.2014', 
          'NC.P.3.2014', 
          'NC.P.5.2014', 
          'NC.P.7.2014', 
          'LC.P.17.2015', 
          'LC.P.23.2015', 
          'NC.P.65.2016', 
          'NC.P.71.2016', 
          'NC.P.72.2016', 
          'NC.P.78.2016', 
          'LC.P.18.2015',
          'CRC.P.22.2015')
          
rrbs.mhl.file <- paste0(data_dir, "/RRBS_180118.gethaplo.mhl.mhbs1.0.useSampleID.txt.gz")
rrbs.umhl.file <- paste0(data_dir, "/RRBS_180118.gethaplo.umhl.mhbs1.0.useSampleID.txt.gz")
rrbs.mhl.dat <- read.table(gzfile(rrbs.mhl.file), sep = '\t', header = T, row.names = 1)
rrbs.umhl.dat <- read.table(gzfile(rrbs.umhl.file), sep = '\t', header = T, row.names = 1)

# add labels to distinguish MHBs with MHL vs uMHL values
mhl_marker_names <- paste0("MHL:",rownames(rrbs.mhl.dat))
umhl_marker_names <- paste0("uMHL:",rownames(rrbs.umhl.dat))

rownames(rrbs.mhl.dat) <- mhl_marker_names
rownames(rrbs.umhl.dat) <- umhl_marker_names
rrbs.combined.dat <- rbind(rrbs.mhl.dat, rrbs.umhl.dat)

rrbs.meta.file <- paste0(data_dir, "/RRBS.getHaplo.sampleInfo.txt")

rrbs.meta <- read.table(rrbs.meta.file, header = T, sep = '\t', stringsAsFactors = F)

rrbs.meta$Cancer.Origin <- rrbs.meta$Tissue ## rename this

rownames(rrbs.meta) <- rrbs.meta$Sample
rrbs.meta <- rrbs.meta[ ! rownames(rrbs.meta) %in% OUTLIERS, ]
intersect.names <- sort(intersect(colnames(rrbs.combined.dat), rownames(rrbs.meta)))
 
rrbs.combined.dat <- rrbs.combined.dat[, intersect.names]
rrbs.meta <- rrbs.meta[intersect.names, ]

print(dim(rrbs.combined.dat))
print(dim(rrbs.meta))

rrbs.plasma.meta <- rrbs.meta[(rrbs.meta$Type == 'Plasma'),]

intersect.names <- sort(intersect(colnames(rrbs.combined.dat), rrbs.plasma.meta$Sample))

#############    
# feature selection & validation for cancer detection in plasma
#############
MHB.filtered.set.file = paste0(result_dir, "/WGBS_MHBv1_MHL_low_in_blood_9-tissue_with_cancer_GSI.list.txt")
MHB.filtered.set <- read.table(MHB.filtered.set.file, sep = '\t', header = T, row.names = 1)

MHB.selected.set <- MHB.filtered.set[MHB.filtered.set$GSI > 0.30,]

print(paste0("Filtered features: ", nrow(MHB.selected.set)))
write.table(MHB.selected.set, file = paste0(result_dir, "/selected.MHBv1.txt"), quote = F, sep = '\t')

rrbs.working.dat <- rrbs.combined.dat[rownames(MHB.selected.set), rrbs.plasma.meta$Sample]
print(dim(rrbs.working.dat))

# Dimension reduction: generate a cluster-level data matrix from the MHB x sample matrix
unique_cluster_list <- unique(MHB.selected.set$group)
length(unique_cluster_list)
rrbs.cluster_level.dat<- matrix(nrow=length(unique_cluster_list), ncol=ncol(rrbs.working.dat))
rrbs.cluster_level.N<- matrix(nrow=length(unique_cluster_list), ncol=ncol(rrbs.working.dat))

for(index in 1:length(unique_cluster_list)){
  all_MHBs_for_one_cluster <- MHB.selected.set[MHB.selected.set$group == unique_cluster_list[index],]
  all_MHBs_for_one_cluster.dat <- rrbs.working.dat[rownames(all_MHBs_for_one_cluster), ] 
  all_MHBs_for_one_cluster.GSI_weighted.dat <- all_MHBs_for_one_cluster.dat*log10(all_MHBs_for_one_cluster$GSI*9)
  rrbs.cluster_level.dat[index,] <- apply(all_MHBs_for_one_cluster.GSI_weighted.dat, 2, function(x) mean(x, na.rm = TRUE))
  rrbs.cluster_level.N[index,] <- apply(all_MHBs_for_one_cluster.dat, 2, function(x) { sum(! is.na(x))})
}

rownames(rrbs.cluster_level.dat) <- unique_cluster_list
colnames(rrbs.cluster_level.dat) <- colnames(rrbs.working.dat)

rrbs.cluster_level.observed.marker <- apply(rrbs.cluster_level.dat, 1, function(x) { sum(! is.na(x))})

rrbs.cluster_level.filtered1.dat <- rrbs.cluster_level.dat[ apply(rrbs.cluster_level.dat, 1, function(x) sum(! is.na(x))/length(x) == 1), ]
rrbs.cluster_level.observed.sites <- apply(rrbs.cluster_level.filtered1.dat, 2, function(x) { sum(! is.na(x))})
#hist(rrbs.cluster_level.observed.sites)

rrbs.cluster_level.selected.dat <- t(rrbs.cluster_level.filtered1.dat)
rrbs.cluster_level.combined.dat <- cbind(rrbs.cluster_level.selected.dat,rrbs.plasma.meta[rownames(rrbs.cluster_level.selected.dat),])

rrbs.cluster_level.combined.dat[,'Type'] <- NULL
rrbs.cluster_level.combined.dat[,'Tissue'] <- NULL
rrbs.cluster_level.combined.dat[,'Sample'] <- NULL
rrbs.cluster_level.combined.dat[,'Cancer.Origin'] <- NULL

trainSet <- rrbs.cluster_level.combined.dat[ rownames(rrbs.cluster_level.combined.dat) %in% normal.vs.cancer.train, ]
testSet <- rrbs.cluster_level.combined.dat[ rownames(rrbs.cluster_level.combined.dat) %in% normal.vs.cancer.test, ]
print(paste0("Normal vs Cancer RF features: ", ncol(trainSet)-1))

set.seed(92976)
metric <- "Accuracy"
ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "grid", classProbs = T, savePredictions = T)
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(Disease ~ . , data = trainSet,  preProc = c("center", "scale"), method = "rf", metric = metric, tuneGrid = tunegrid, trControl = ctrl)

#validate on held-out samples
result.predicted <- predict(rf_gridsearch, testSet) # Prediction
confusionMatrix(result.predicted, testSet$Disease)

result.predicted.prob <- predict(rf_gridsearch, testSet, type='prob') # Prediction
roc.res <- roc(ifelse(testSet$Disease=="Cancer",1,0), result.predicted.prob[[2]])
print(paste0("AUC=", roc.res$auc))
pdf(paste0(result_dir, "/Normal_vs_cancer_75pct_split_AUC_plot.pdf"), width=3, height=3)
plot.roc(roc.res, main = "Normal vs Cancer")
dev.off()

for(seed_index in 1:length(seed.table)){
  set.seed(seed.table[seed_index])
  index <- createDataPartition(rrbs.cluster_level.combined.dat$Disease, p=0.75, list=FALSE)

  trainSet <- rrbs.cluster_level.combined.dat[ index,]
  testSet <- rrbs.cluster_level.combined.dat[-index,]
  print(paste0("Normal vs Cancer RF features: ", ncol(trainSet)-1))

  metric <- "Accuracy"
  ctrl <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid", classProbs=T, savePredictions = T)
  tunegrid <- expand.grid(.mtry=c(1:15))
  rf_gridsearch <- train(Disease ~ . , data = trainSet,  preProc = c("center", "scale"), method="rf", metric=metric, tuneGrid=tunegrid, trControl=ctrl)

  #validate on held-out samples
  result.predicted <- predict(rf_gridsearch, testSet) # Prediction
  confusionMatrix(result.predicted, testSet$Disease)    

  result.predicted.prob <- predict(rf_gridsearch, testSet, type='prob') # Prediction
  roc.res <- roc(ifelse(testSet$Disease=="Cancer",1,0), result.predicted.prob[[2]])
  print(paste0("AUC=", roc.res$auc))
}


#############
# feature selection & validation for separating lung versus colon cancer plasma
#############


MHB.filtered.set.file = paste0(result_dir, "/WGBS_MHBv1_MHL_low_in_blood_9-tissue_with_cancer_GSI.list.txt")
MHB.filtered.set <- read.table(MHB.filtered.set.file, sep = '\t', header = T, row.names = 1)

MHB.selected.set <- MHB.filtered.set[MHB.filtered.set$GSI > 0.45,]
MHB.selected.set <- MHB.selected.set[ - grep("cancer", MHB.filtered.set$group),] # remove all features related to cancers

print(paste0("Filtered features: ", nrow(MHB.selected.set)))
write.table(MHB.selected.set, file=paste0(result_dir, "/selected.MHBv1.txt"), quote = F, sep = '\t')

rrbs.cancer.lung.meta <- rrbs.plasma.meta[(rrbs.plasma.meta$Disease == 'Cancer')&(rrbs.plasma.meta$Cancer.Origin == 'Lung'),]
rrbs.cancer.colon.meta <- rrbs.plasma.meta[(rrbs.plasma.meta$Disease == 'Cancer')&(rrbs.plasma.meta$Cancer.Origin == 'Colon'),]
rrbs.plasma.selected.meta <- rbind(rrbs.cancer.lung.meta, rrbs.cancer.colon.meta)

rrbs.working.dat <- rrbs.combined.dat[rownames(MHB.selected.set), rrbs.plasma.selected.meta$Sample]
print(dim(rrbs.working.dat))

# Dimension reduction: generate a cluster-level data matrix from the MHB x sample matrix
unique_cluster_list <- unique(MHB.selected.set$group)
length(unique_cluster_list)
rrbs.cluster_level.dat <- matrix(nrow=length(unique_cluster_list), ncol=ncol(rrbs.working.dat))
rrbs.cluster_level.N <- matrix(nrow=length(unique_cluster_list), ncol=ncol(rrbs.working.dat))

for(index in 1:length(unique_cluster_list)){
  all_MHBs_for_one_cluster <- MHB.selected.set[MHB.selected.set$group == unique_cluster_list[index],]
  all_MHBs_for_one_cluster.dat <- rrbs.working.dat[rownames(all_MHBs_for_one_cluster), ] 
  all_MHBs_for_one_cluster.GSI_weighted.dat <- all_MHBs_for_one_cluster.dat*log10(all_MHBs_for_one_cluster$GSI*9)
  rrbs.cluster_level.dat[index,] <- apply(all_MHBs_for_one_cluster.GSI_weighted.dat, 2, function(x) mean(x, na.rm = TRUE))
  rrbs.cluster_level.N[index,] <- apply(all_MHBs_for_one_cluster.dat, 2, function(x) { sum(! is.na(x))})
}

rownames(rrbs.cluster_level.dat) <- unique_cluster_list
colnames(rrbs.cluster_level.dat) <- colnames(rrbs.working.dat)

rrbs.cluster_level.observed.marker <- apply(rrbs.cluster_level.dat, 1, function(x) { sum(! is.na(x))})

rrbs.cluster_level.filtered1.dat <- rrbs.cluster_level.dat[ apply(rrbs.cluster_level.dat, 1, function(x) sum(! is.na(x))/length(x) == 1), ]
rrbs.cluster_level.observed.sites <- apply(rrbs.cluster_level.filtered1.dat, 2, function(x) { sum(! is.na(x))})

rrbs.cluster_level.selected.dat <- t(rrbs.cluster_level.filtered1.dat)
rrbs.cluster_level.combined.dat <- cbind(rrbs.cluster_level.selected.dat,rrbs.plasma.meta[rownames(rrbs.cluster_level.selected.dat),])

rrbs.cluster_level.combined.dat[,'Type'] <- NULL
rrbs.cluster_level.combined.dat[,'Disease'] <- NULL
rrbs.cluster_level.combined.dat[,'Sample'] <- NULL
rrbs.cluster_level.combined.dat[,'Cancer.Origin'] <- NULL

trainSet <- rrbs.cluster_level.combined.dat[ rownames(rrbs.cluster_level.combined.dat) %in% colon.vs.lung.train, ]
testSet <- rrbs.cluster_level.combined.dat[ rownames(rrbs.cluster_level.combined.dat) %in% colon.vs.lung.test, ]
print(paste0("Colon vs Lung RF features: ", ncol(trainSet)-1))

set.seed(252241)
metric <- "Accuracy"
ctrl <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid", classProbs = T, savePredictions = T)
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(Tissue ~ . , data = trainSet,  preProc = c("center", "scale"), method = "rf", metric = metric, tuneGrid = tunegrid, trControl = ctrl)

#validate on held-out samples
result.predicted <- predict(rf_gridsearch, testSet) # Prediction
confusionMatrix(result.predicted, testSet$Tissue)

result.predicted.prob <- predict(rf_gridsearch, testSet, type = 'prob') # Prediction
roc.res <- roc(ifelse(testSet$Tissue == "Lung", 1, 0), result.predicted.prob[[2]])
print(paste0("AUC=", roc.res$auc))
pdf(paste0(result_dir, "/Colon_vs_lung_75pct_split_AUC_plot.pdf"), width = 3, height = 3)
plot.roc(roc.res, main = "Colon vs Lung")
dev.off()


for(seed_index in 1:length(seed.table)){
  set.seed(seed.table[seed_index])
  index <- createDataPartition(rrbs.cluster_level.combined.dat$Tissue, p = 0.75, list = FALSE)
   
  trainSet <- rrbs.cluster_level.combined.dat[ index,] 
  testSet <- rrbs.cluster_level.combined.dat[ - index,] 
  print(paste0("Colon vs Lung RF features: ", ncol(trainSet) - 1))
    
  metric <- "Accuracy"
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "grid", classProbs  = T, savePredictions = T)
  tunegrid <- expand.grid(.mtry=c(1:15))
  rf_gridsearch <- train(Tissue ~ . , data = trainSet,  preProc = c("center", "scale"), method = "rf", metric = metric, tuneGrid = tunegrid, trControl = ctrl)

  #validate on held-out samples
  result.predicted <- predict(rf_gridsearch, testSet) # Prediction
  print(confusionMatrix(result.predicted, testSet$Tissue))    

  result.predicted.prob <- predict(rf_gridsearch, testSet, type='prob') # Prediction
  roc.res <- roc(ifelse(testSet$Tissue == "Lung", 1, 0), result.predicted.prob[[2]])
  print(paste0("AUC=", roc.res$auc))
}

