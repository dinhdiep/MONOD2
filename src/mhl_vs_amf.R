library(ggplot2)
library(compiler)

#### GSI function ####
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

## the function can be be viewed as a two step process
## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
## using simple options or by a user supplied function
## 2. with the now resahped data the plot, the chosen labels and plot style are built
ggheat <- function(m, rescaling = 'none', clustering = 'none', labCol = T, labRow = T, border = F, 
                   heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato'), legend.title = NULL, title="Heatmap") {
  require(reshape)
  require(ggplot2)
  
  ## you can either scale by row or column not both! 
  ## if you wish to scale by both or use a differen scale method then simply supply a scale
  ## function instead NB scale is a base funct
  
  if(is.function(rescaling))
  { 
    m=rescaling(m)
  } 
  else 
  {
    if(rescaling=='column') 
      m=scale(m, center=T)
    if(rescaling=='row') 
      m=t(scale(t(m),center=T))
  }
  
  ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function 
  
  if(is.function(clustering)) 
  {
    m=clustering(m)
  }else
  {
    if(clustering=='row')
      m=m[hclust(dist(m))$order, ]
    if(clustering=='column')  
      m=m[,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
  }
  ## this is just reshaping into a ggplot format matrix and making a ggplot layer
  
  rows=dim(m)[1]
  cols=dim(m)[2]
  melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(m))
  g=ggplot(data=melt.m)
  
  ## add the heat tiles with or without a white border for clarity
  
  if(border==TRUE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='white')
  if(border==FALSE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
  
  ## add axis labels either supplied or from the colnames rownames of the matrix
  
  if(labCol==T) 
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m), expand = c(0.005,0))
  if(labCol==F) 
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
  if(labRow==T)
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m), expand = c(0.005,0))	
  if(labRow==F) 
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
  
  ## get rid of grey panel background and gridlines
  
  g2 = g2 + theme(panel.grid.minor = element_line(colour=NA), panel.grid.major = element_line(colour=NA),
                  panel.background = element_rect(fill=NA, colour=NA), axis.text.x = element_text(angle = 90, hjust = 1),
                  axis.ticks = element_blank()) + ggtitle(title)
  
  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  return(g2+scale_fill_gradient2(low = heatscale[1], mid = heatscale[2], high = heatscale[3], guide = guide_colorbar(title = legend.title)))
}


# Load datasets
wgbs.mhl.file <- "data/ng.3805/WGBS.adult_tissues.mhl.txt.gz"
wgbs.amf.file <- "data/ng.3805/WGBS.adult_tissues.amf.txt.gz"
dmr.bed.file <- "data/ng.3805/RRBS_MHBs.sorted.DMR.withID.bed"

wgbs_mhl_data <- read.table(gzfile(wgbs.mhl.file), sep = "\t", header = T, row.names = 1)
wgbs_amf_data <- read.table(gzfile(wgbs.amf.file), sep = "\t", header = T, row.names = 1)

pdf("plots.pdf", width = 5, height = 5)

# Define samples we want to use
adult_stl_samples <- c('STL002SB.01', 'STL003SB.01', 
                       'STL003SG.01', 'STL002EG.01', 
                       'STL003EG.01', 'STL002FT.01', 
                       'STL003FT.01', 'STL003LV.01', 
                       'STL003RV.01', 'STL002AD.01', 
                       'STL003AD.01', 'STL011LI.01', 
                       'STL002LG.01', 'STL002PO.01', 
                       'STL003PO.01', 'STL002OV.01', 
                       'STL002PA.01', 'STL003PA.01', 
                       'STL002SX.01', 'STL003SX.01', 
                       'STL002GA.01', 'STL003GA.01', 
                       'STL002AO.01', 'STL003AO.01', 
                       'STL003RA.01')

# Define sample types
adult_stl_labels <- c('colon 1', 'colon 2', 'colon 3', 'esophagus 1', 'esophagus 2', 'fat 1', 'fat 2', 'heart 1', 
                      'heart 2', 'kidney 1', 'kidney 2', 'liver 1', 'lung 2', 'muscle 1', 'muscle 2', 'ovary 1', 
                      'pancreas 1', 'pancreas 2', 'spleen 1', 'spleen 2', 'stomach 1', 'stomach 2', 'vessel 1', 
                      'vessel 2', 'vessel 3')

adult_stl_labels <- sapply(adult_stl_labels, function(x) strsplit(x, split = " ")[[1]][[1]])

# Load DMR MHB regions
dmr_mhbs <- read.table(dmr.bed.file, sep = "\t", header = F, stringsAsFactors = F)
colnames(dmr_mhbs) <- c("chrom", "start", "end", "Region")
dmr_mhb_list <- dmr_mhbs$Region

# Subset MHL matrix to DMR regions and STL samples
wgbs_mhl_data <- wgbs_mhl_data[dmr_mhb_list, adult_stl_samples]
colnames(wgbs_mhl_data) <- adult_stl_labels

# Subset AMF matrix to DMR regions and STL samples
wgbs_amf_data <- wgbs_amf_data[dmr_mhb_list, adult_stl_samples]
colnames(wgbs_amf_data) <- adult_stl_labels

wgbs_mhl_data <- wgbs_mhl_data[apply(wgbs_mhl_data, 1, function(x) sum(! is.na(x))/length(x) >= 0.8),]
wgbs_amf_data <- wgbs_amf_data[rownames(wgbs_mhl_data),]

# Calculate GSI using the AMF matrix
gsi.df <- gsi(wgbs_amf_data, min.frac = 1)

# Take top 150 DMRs for each tissue
gsi.df <- gsi.df[order(gsi.df$GSI, decreasing = T),]
top.gsi.df <- Reduce(rbind, by(gsi.df, gsi.df["group"], head, n = 150))

# Generate heatmaps
ggheat(as.matrix(wgbs_mhl_data[top.gsi.df$region,]), labCol = T, labRow = F, title="MHL")
ggheat(as.matrix(wgbs_amf_data[top.gsi.df$region,]), labCol = T, labRow = F, title="AMF")

# Quantitative signal to noise analysis by looking at ratio of diagonal to off-diagonal values
amf_snr <- rep(0, nrow(top.gsi.df))
mhl_snr <- rep(0, nrow(top.gsi.df))

wgbs_mhl_data <- as.matrix(wgbs_mhl_data)
wgbs_amf_data <- as.matrix(wgbs_amf_data)

# Calculate signal to noise for each region
for (i in 1:nrow(top.gsi.df)) {
  region <- top.gsi.df[i, "region"]
  group <- top.gsi.df[i, "group"]
  
  mhl.fg <- mean(wgbs_mhl_data[region, which(colnames(wgbs_mhl_data) == group)], na.rm = T)
  mhl.bg <- mean(wgbs_mhl_data[region, which(colnames(wgbs_mhl_data) != group)], na.rm = T)
  
  amf.fg <- mean(wgbs_amf_data[region, which(colnames(wgbs_amf_data) == group)], na.rm = T)
  amf.bg <- mean(wgbs_amf_data[region, which(colnames(wgbs_amf_data) != group)], na.rm = T)
  
  amf_snr[[i]] <- amf.fg/amf.bg
  mhl_snr[[i]] <- mhl.fg/mhl.bg
}

# Visualize MHL signal to noise vs AMF signal to noise
gg.df <- data.frame(mhl_snr, amf_snr)
ggplot(gg.df, aes(amf_snr, mhl_snr)) + geom_point(col='slateblue') +
  theme_classic() + theme(legend.position="none", text = element_text(size = 12)) +
  xlab("AMF Signal to Noise Ratio") + ylab("MHL Signal to Noise Ratio") + geom_abline(slope = 1, intercept = 0, col='seagreen4') + 
   scale_x_log10(
   lim = c(0.1, 1.1e3),
   breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
 scale_y_log10(
   lim = c(0.1, 1.1e3),
   breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x)))

dev.off()
