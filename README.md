# MONOD2
MONOD2 is a toolkit for methylation haplotype analysis for bisulfite sequencing or (BS-seq) data using high-throughput sequencing.  The analysis can be divided into three parts. First sequence alignment files are analyzed to generate methylation haplotypes which includes haplotype strings and haplotype counts for CpG positions with coverage.  Next, a collection of methylation haplotype files can be analyzed for identification of methylation haplotype blocks or for differential methylation haplotype load. 

MONOD2 is comprised of executable 'sh', 'perl', and 'R' programs.

## Contact
email: hdinhdp@gmail.com

## Download

You can use git to download the entire toolkit. 

'''
Code
'''

## Required R packages

'''
require(ggplot2)
require(earth)
require(skmeans)
require(ModelMetrics)
require(caret)
require(ROCR)
require(gplot2)
require(reshape)
require(compiler)

'''

## Required pre-installed software

1. samtools version 1.2 or above
2. bedtools version 2.26 or above 

## Important Notes

MONOD2 is based on the approach described in [! Guo et al. 2017 Nature Genetics][doi:10.1038/ng.3805]. Here, we have made significant modifications from the original methods and codes that was hosted on the [! Supplementary Website][http://genome-tech.ucsd.edu/public/MONOD_NG_TR44413]. 
