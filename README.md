# MONOD2
`MONOD2` is a toolkit for methylation haplotype analysis for bisulfite sequencing or (BS-seq) data using high-throughput sequencing.  The analysis can be divided into three parts. First sequence alignment files are analyzed to generate methylation haplotypes which includes haplotype strings and haplotype counts for CpG positions with coverage.  Next, a collection of methylation haplotype files can be analyzed for identification of methylation haplotype blocks or for differential methylation haplotype load. 

MONOD2 is comprised of executable `sh`, `perl`, and `R` programs. No installation is necessary except for the required dependencies listed below.

## Contact
email: hdinhdp@gmail.com

## Download

You can use git to download the entire toolkit. 

```
git clone https://github.com/dinhdiep/MONOD2

```

## Required R packages

You can install all the R packages from command line in R with the following:

```
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

## Required pre-installed software

1. [samtools version 1.2 or above](http://samtools.sourceforge.net/)
2. [bedtools version 2.26 or above](http://bedtools.readthedocs.io/en/latest/)

## Important Notes

`MONOD2` is based on the idea described in Guo et al. 2017 Nature Genetics (doi:10.1038/ng.3805). Here, we have made significant modifications from the original methods and codes that was hosted on the [Supplementary Website](http://genome-tech.ucsd.edu/public/MONOD_NG_TR44413). Some notable modifications are listed below:

1. The original code for Figure 2 was not producing the correct values from the Figure. Note that there are two errors in Figure 2 as follows: the epipolymorphism for panel 4 should be 0.9375 (the published figure have 0.375 due to mis-editing), the MHL value for panel 5 should be 0.1167 (the published figure have a rounded up value of 0.1200).  

2. The original comparison of AMF versus MHL at tissue specific regions were unfair because they used differentially methylated regions identified by MHL. In `MONOD2` we use only the MHB regions overlapping with published tissue specific DMRs to compare AMF versus MHL.

3. The original method for tumor load estimation performed features selection without separation of test data (the plasma samples). In `MONOD2` we identified the features using a subset of 50 normal plasma only, and test on the remaining normal and cancer patient plasma samples.

4. The original method for methylation haplotype load analysis of cf-DNA did not have proper separation of training and test data for features selection although in model building the training and test data were separated. We realized that this would heavily bias the results to the sample set. In `MONOD2` we separated training and test data for both features selection and model building.

## Usage

### Extracting methylation haplotypes from sequence alignment files

Sequence alignment files in BAM format should have the expected sam flags for Watson (forward) and Crick (reverse) strands. The code below will generate CpG haplotype files without clonal removal (any clonal remove must be done upstream).
```
bam2cghap.sh [cpg position file] [bam file] [output file prefix name]

```

Run example
```
./scripts/bam2cghap.sh allcpg/cpg.small.txt.gz BAMfiles/Colon_primary_tumor_sept9_promoter.bam test

```
The version 1 code is provided but user must specify RRBS or WGBS mode and can only use BAM that have been mapped to hg19. For RRBS and WGBS data are treated differently with the major difference in that reads which are considered clonal in WGBS would be removed in order to report only one haplotype.

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

Methylation haplotypes were split into continuous mappable bins and then an algorithm greedily select the largest possible continuous region with a minimum linkage disequilibrium score to be considered methylation haplotype blocks. Blocks must have at least 3 CpGs sites. For more info about methylation linkage disequilibrium, please see (Shoemaker et al Genome Research)[http://genome.cshlp.org/content/20/7/883.full.pdf].

```
cghap2mhbs.sh [haplotype file] [target bed] [minimum LD R2 cutoff] [output name prefix]

```
Run example

```
./scripts/cghap2mhbs.sh HaploInfo/chr22.sub.hapInfo.txt example/N37_10_tissue_pooled.autosomes.RD10_80up.genomecov.bed 0.3 chr22

```

### Compare AMF versus MHL 

We can compare the signal to noise levels between a matrix calculated using the average methylation frequency (AMF) and a matrix calculated using the methylation haplotype load (MHL). The three required files are as followed.
1. WGBS MHL matrix. Example: ng.3805/WGBS.getHaplo.mhl.mhbs1.0.rmdup_consistent.useSampleID.txt.gz
2. WGBS AMF matrix. Example: ng.3805/WGBS.getHaplo.amf.mhbs1.0.rmdup_consistent.useSampleID.txt.gz
3. The published DMRs file. Example: ng.3805/RRBS_MHBs.sorted.DMR.withID.bed

The only requirements for (1) and (2) are that the region indices and the sample IDs are the same between the two matrices. The region file (3) is a list of regions (subset of the region indices used to make the matrices). 

Run example

```
./scripts/mhl_vs_amf.sh

```

The output is an Rplots.pdf file.

### Tumor load estimation

First, we generate the simulation files by merging a subset of normal plasma data with cancer tissue data. We have both primary tumor tissue data for lung cancer and colon cancer. Since the large files cannot be hosted on github, we randomly sampled (with replacement) the datasets and provided them in the BAMfiles folder. Note that read IDs can contain either mapped reads or all reads. We performed the analysis on all reads so that the mappability of sampled fragments also reflect the mappability of 'non-cancer' versus 'cancer' fragments. 

The script needs to be modified from within (using a text editor) in order to run properly. The following lines must be edited.

```
numsimulation=20 # for testing use 2
n1=1000000 # for testing use 10000

cctfastq="BAMfiles/CCT.readIDs.txt"
lctfastq="BAMfiles/LCT.readIDs.txt"
ncpfastq="BAMfiles/NCP.readIDs.txt"

cctbam="BAMfiles/CCT.bam"
lctbam="BAMfiles/LCT.bam"
ncpbam="BAMfiles/NCP.bam"

allcpgs="allcpg/hg19.fa.allcpgs.txt.gz"

```
Run example

```
 ./scripts/simulate-cghap.sh

```

After all the simulation files have been generated, we can run features selection and tumor load estimation using the simulated data files.

### Plasma prediction



