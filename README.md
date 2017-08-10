# MONOD2
MONOD2 is a toolkit for methylation haplotype analysis for bisulfite sequencing or (BS-seq) data using high-throughput sequencing.  The analysis can be divided into three parts. First sequence alignment files are analyzed to generate methylation haplotypes which includes haplotype strings and haplotype counts for CpG positions with coverage.  Next, a collection of methylation haplotype files can be analyzed for identification of methylation haplotype blocks or for differential methylation haplotype load. Methylation haplotype blocks are regions with high linkage disequilibrium between pairs of CpGs and they were found to be enriched at transcription factor binding sites and enhancers. For differentialy methylation haplotype load analysis we provide two examples where tumor load estimation and plasma tissue of origin prediction can be performed on cf-DNA.  

MONOD2 is comprised of executable `bash`, `perl`, and `R` programs. No installation is necessary except for the required dependencies listed below.

## Contact
email: hdinhdp@gmail.com

## Download

You can use git to download the entire codebase and datasets. 

```
git clone https://github.com/dinhdiep/MONOD2

```

## Required R packages

You can install all the R packages from command line in R with the following code.

```
require(gplots)
require(MASS)
require(pheatmap)
require(RColorBrewer)
require(randomForest)
require(impute)
require(abind)
require(ggplot2)
require(earth)
require(skmeans)
require(ModelMetrics)
require(caret)
require(ROCR)
require(gplot2)
require(reshape)
require(compiler)

```

## Required software

The following two software must be installed.

1. [samtools version 1.2 or above](http://samtools.sourceforge.net/)
2. [bedtools version 2.26 or above](http://bedtools.readthedocs.io/en/latest/)

## Important Notes

MONOD2 is based on the idea described in Guo et al. 2017 Nature Genetics (doi:10.1038/ng.3805). However, we have made modifications from the original methods and codes that was originally hosted on the [Supplementary Website](http://genome-tech.ucsd.edu/public/MONOD_NG_TR44413). A copy of the original processed data, code, and supplementary tables are available here as well in the ng.3805 directory. Some notable modifications are listed below:

1. The original code for Figure 2 was not producing the correct values from the Figure. Note that there are two errors in Figure 2 as follows: the epipolymorphism for panel 4 should be 0.9375 (the published figure have 0.375 due to mis-editing), the MHL value for panel 5 should be 0.1167 (the published figure have a rounded up value of 0.1200).  

2. The original comparison of AMF versus MHL at tissue specific regions were unfair because they used differentially methylated regions identified by MHL. In MONOD2 we use the MHB regions overlapping with published tissue specific DMRs to make the comparison.

3. The original method for tumor load estimation performed features selection without separation of test data (the plasma samples). In MONOD2 we identified the features using a subset of 50 normal plasma, and then computed tumor load on the remaining normal and cancer patient plasma samples.

4. The original method for classification of cf-DNA did not have proper separation of training and test data during features selection although in model building the training and test data were separated. We realized that this would heavily bias the results to the sample set. In MONOD2 we separated training and test data for both features selection and model building.

## Usage

### Extracting methylation haplotypes from sequence alignment files

Sequence alignment files in BAM format should have the expected sam flags for Watson (forward) and Crick (reverse) strands. The code below will generate CpG haplotype files without clonal removal (any clonal removal must be performed upstream).

```
bam2cghap.sh [cpg position file] [bam file] [output file prefix name]

```

Run example

```
./scripts/bam2cghap.sh allcpg/cpg.small.txt.gz BAMfiles/Colon_primary_tumor_sept9_promoter.bam test

```
Version 1 code is provided but user must specify RRBS or WGBS mode and can only use BAM that have been mapped to hg19. RRBS and WGBS data are treated differently in that reads which are considered clonal in WGBS would be removed.

Run RRBS example

```
./scripts/bam2cghap_v1.sh RRBS Colon_primary_tumor_sept9_promoter.RD1_80up.genomecov.bed allcpg/cpg.small.txt.gz BAMfiles/Colon_primary_tumor_sept9_promoter.bam test

```
Run WGBS example

```
./scripts/bam2cghap_v1.sh WGBS Colon_primary_tumor_sept9_promoter.RD1_80up.genomecov.bed allcpg/cpg.small.txt.gz BAMfiles/Colon_primary_tumor_sept9_promoter.bam test

```

### Making mappable bin file

Empirically determine which regions in the genome have high mappability using whole genome bisulfite sequencing datasets. 

```
make-mappable-bins.sh [bam file] [minimum depth cutoff]

```

Run example

```
./scripts/make-mappable-bins.sh BAMfiles/Colon_primary_tumor_sept9_promoter.bam 5

```

### Identifying the methylation haplotype blocks 

Methylation haplotypes are split into continuous mappable bins and then an algorithm greedily selects the largest possible continuous region with a minimum linkage disequilibrium score to be considered methylation haplotype blocks. Blocks must have at least 3 CpGs sites. The original idea for methylation haplotype blocks were described in details by [Shoemaker et al. 2010 Genome Research](http://genome.cshlp.org/content/20/7/883.full.pdf).

```
cghap2mhbs.sh [haplotype file] [target bed] [minimum LD R2 cutoff] [output name prefix]

```
Run example

```
./scripts/cghap2mhbs.sh HaploInfo/chr22.sub.hapInfo.txt example/N37_10_tissue_pooled.autosomes.RD10_80up.genomecov.bed 0.3 chr22

```

### Compare AMF versus MHL 

We can compare the signal to noise levels between a matrix calculated using the average methylation frequency (AMF) and a matrix calculated using the methylation haplotype load (MHL) to demonstrate the advantage of using MHL on heterogenous samples. The following files are required.

1. WGBS MHL matrix. Example: ng.3805/WGBS.getHaplo.mhl.mhbs1.0.rmdup_consistent.useSampleID.txt.gz
2. WGBS AMF matrix. Example: ng.3805/WGBS.getHaplo.amf.mhbs1.0.rmdup_consistent.useSampleID.txt.gz
3. List of MHBs overlapping with published DMRs. Example: ng.3805/RRBS_MHBs.sorted.DMR.withID.bed

The only requirements for (1) and (2) are that the region indices and the sample IDs are the same between the two matrices. The region file (3) is a list of regions (subset of the region indices used to make the matrices). 

Run example

```
Rscript src/Rcode/mhl_vs_amf.R ng.3805/WGBS.getHaplo.mhl.mhbs1.0.rmdup_consistent.useSampleID.txt.gz ng.3805/WGBS.getHaplo.amf.mhbs.rmdup_consistent.useSampleID.txt.gz ng.3805/RRBS_MHBs.sorted.DMR.withID.bed

```

The output is an Rplots.pdf file that includes two heatmaps similar to Figure 3 from Guo et al. 2017, while subsetting at regions with tissue specific DMRs and also includes a scatteplot showing the relative signal to noise levels for AMF versus MHL.


### Preprocess the data matrices 

Analysis cannot be performed across RRBS and WGBS datasets due to techical artifacts that may bias the results. Therefore we generate two separate matrices, one for each dataset. We then perform data pruning where samples with too few region coverage and regions with too few sample coverage are removed. After pruning, imputation using k-nearest neighbor to fill all the missing values is performed. To run the provided preprocessing script, the following files are required.

1. The metadata for the WGBS tissue samples: ng.3805/WGBS.getHaplo.sampleInfo.txt
2. The metadata for the RRBS tissue/plasma samples: ng.3805/RRBS.getHaplo.sampleInfo.txt
3. WGBS MHL matrix: ng.3805/WGBS.getHaplo.mhl.mhbs1.0.rmdup_consistent.useSampleID.txt
4. RRBS MHL matrix: ng.3805/RRBS_170609.gethaplo.mhl.mhbs1.0.useSampleID.txt

Usage info

```
preprocess.sh [rrbs.matrix.gz] [wgbs.matrix.gz] [rrbs.meta] [wgbs.meta] [output Rdata name]

```

Run example

```
./scripts/preprocess.sh ng.3805/RRBS_170609.gethaplo.mhl.mhbs1.0.useSampleID.txt.gz ng.3805/WGBS.getHaplo.mhl.mhbs1.0.rmdup_consistent.useSampleID.txt.gz ng.3805/RRBS.getHaplo.sampleInfo.txt ng.3805/WGBS.getHaplo.sampleInfo.txt wgbs_rrbs_clean.Rdata

```

The output is an R data file named 'wgbs_rrbs_clean.Rdata' and several PDF files showing the data quality.

### Tumor load estimation

First, we generate the simulation files by merging a subset of normal plasma data with cancer tissue data. We have both primary tumor tissue data for lung cancer and colon cancer. Since the large files cannot be hosted on github, we randomly sampled (with replacement) the datasets and provided them in the BAMfiles folder. Note that read IDs can contain either mapped reads or all reads. We performed the analysis on all reads so that the mappability of sampled fragments also reflect the mappability of 'non-cancer' versus 'cancer' fragments. 

The script needs to be modified from within (using a text editor) in order to run user-provided read IDs and BAM files. The following lines as is will perform the analysis on the example read IDs and BAM files. 

```
numsimulation=20 # number of simulations to perform
n1=1000000 # number of reads per simulation 
 
cctfastq="BAMfiles/CCT.readIDs.txt" # path to the read IDs for colon tumor samples
lctfastq="BAMfiles/LCT.readIDs.txt" # path to the read IDs for lung tumor samples
ncpfastq="BAMfiles/NCP.readIDs.txt" # path to the read IDs for normal plasma training subset

cctbam="BAMfiles/CCT.bam" # path to the merged BAM for colon tumor samples
lctbam="BAMfiles/LCT.bam" # path to the merged BAM for lung tumor samples
ncpbam="BAMfiles/NCP.bam" # path to the merged BAM for normal plasma training subset

cpgs="allcpg/hg19.fa.allcpgs.txt.gz" # path to the CpG positions file

```
Run example

```
 ./scripts/simulate-cghap.sh

```

There should be many simulated files generated including two list files.

1. cct.hapinfo.list_20 
2. lct.hapinfo.list_20 

The should also be three directories which contain simulated data for each sample type.
1. CCT_simulation 
2. LCT_simulation 
3. NCP_simulation

Next, we can run features selection and tumor load estimation using the simulated data files. The script for tumor load estimation have the following usage.

```
tumor_load_estimation.sh [rdata file] [colon cancer simulation table] [lung cancer simulation table] [list of normal plasma training samples] [number of simulations]

```

Run example

```
./scripts/tumor_load_estimation.sh wgbs_rrbs_clean.Rdata cct.hapinfo.list_20 lct.hapinfo.list_20 ng.3805/subset_normal_plasma_training 20

```

The output files are as follows. 

1. The Rdata which saves all the outputs : Tumor_load_estimation.Rdata
2. Boxplots similar to Figure 4d from Guo et al. 2017 which shows the estimated tumor proportions for the plasma samples: estimated_proportions.pdf
3. Boxplots similar to Figure 4a,b from Guo et al. 2017 which shows the differential MHL levels in tumor marker regions for different sample types: mhl_boxplot.pdf
4. Heatmaps similar to Figure 4a,b from Guo et al. 2017 which shows the MHL levels in normal plasma, cancer plasma, and cancer tissues in the tumor marker regions: heatmap.pdf
5. A table of the standard curve values generated from simulated data: standard_curves_values.txt
6. A list of lung cancer markers used in the analysis: lung_cancer_markers.txt
7. A list of colon cancer markers used in the analysis: colon_cancer_markers.txt


### Plasma prediction

To perform plasma prediction, we asked if plasma samples from healthy individuals, lung cancer patients, and colon cancer patients can be identified from their methylation haplotype load profiles. In the provided code, a random sample of 20 colon cancer plasma, 20 lung cancer plasma, and 30 healthy plasma were used to train an ensemble MARS (Multivariate Adaptive Regression Splines) model (R `earth` package). Multiclass prediction is performed by a linear discriminant model on prediction scores. The linear discriminant model was built using the prediction scores of training data. 

Run example

```
Rscript src/Rcode/plasma_prediction.R wgbs_rrbs_clean.Rdata

```

The resulting AUCs and confusion matrix are printed to screen while AUC curves are generated in an Rplots.pdf file.

```
$mauc
[1] 0.7770558

$auc
    Colon      Lung    Normal
0.7250000 0.7505669 0.8556005

          Reference
Prediction Colon Lung Normal
    Colon      4    4      2
    Lung       3    6      0
    Normal     4    5     30

```

The model related files

1. List of features. The final ensemble model have 154 unique features only: final.marker.list.txt 
2. Saved Rdata of the final model: final.ensemble.model.Rdata



 
