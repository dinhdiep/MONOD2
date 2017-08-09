#!/usr/bin/R
# Methylation haplotypes to Methylation haplotype load (mode weight i)
# Contact: Shicheng Guo (shihcheng.guo@gmail.com)
# Version 1.4
# Update: May 5 2017

# R package "entropy" is required for calculating methylation entropy. This package implements various estimators of entropy, such as the shrinkage estimator by Hausser and Strimmer, the maximum likelihood and the Millow-Madow estimator, various Bayesian estimators, and the Chao-Shen estimator. It also offers an R interface to the NSB estimator. Furthermore, it provides functions for estimating Kullback-Leibler divergence, chi-squared, mutual information, and chi-squared statistic of independence. In addition there are functions for discretizing continuous random variables. More information at: https://CRAN.R-project.org/package=entropy.


# This software is Copyright � 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu

# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#  
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.

if (!require("entropy")) install.packages("entropy")

hap_comb_tmp<-function(hapinfo.single){
  # hap_com.single is a character of haplotype, such as "TCTCTT"
  hap_com.single<-c()
  for(i in 1:nchar(hapinfo.single)){
    for(j in i:(nchar(hapinfo.single))){
      hap_com.single<-c(hap_com.single,(substr(hapinfo.single,i,j)))
    } 
  }
  return(hap_com.single)
}

hap_comb<-function(hapinfo){
  # haplotype is array of observed methylation haplotype and return
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

mhl<-function(hapinfo){
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb(hapinfo))
  weight=max(unlist(lapply(hapinfo,function(x) length(unlist(strsplit(x,""))))))
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  for(i in 1:weight){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep(paste(rep("C",i),collapse = ""),hap_comb_array[input_tmp]))
    mhl_value<-c(mhl_value,i*nmeth/length(input_tmp))
  }
  mhl_value<-sum(mhl_value)/sum(1:weight)
  mhl_value
  return(mhl_value)
}

mf<-function(hapinfo){
  mf_value<-c()
  hap_comb_array<-unlist(hapinfo)
  meth_hap=unlist(lapply(hap_comb_array,function(x) unlist(strsplit(x,""))))
  nmeth<-length(grep("C",meth_hap))
  nunmeth<-length(grep("T",meth_hap))
  mf_value<-nmeth/(nunmeth+nmeth)
  return(mf_value)
}

methentropy<-function(hapinfo){
  library("entropy")
  entropy<-entropy(table(hapinfo),unit="log2")/nchar(hapinfo[1])
  return(entropy)
}

ppoly<-function(hapinfo){
  ppoly<-1-sum((table(hapinfo)/length(hapinfo))^2)
  return(ppoly)
}

hapinfo1<-c(rep("CCCC",16))
hapinfo2<-c(rep("TTTT",16))
hapinfo3<-c(rep("TTTT",8), rep("CCCC",8))
hapinfo5<-c(rep("CCTT",8), rep("TTCC",8))
hapinfo4<-c("TTTT","TTTC", "TTCT", "TCTT", "CTTT", "TTCC","TCTC","CTTC","TCCT","CTCT","CCTT","TCCC","CTCC","CCTC","CCCT","CCCC")
results<-data.frame(Fig2A=c(100*mf(hapinfo1), methentropy(hapinfo1),ppoly(hapinfo1),length(unique(hapinfo1)), mhl(hapinfo1)),
			Fig2B=c(100*mf(hapinfo2), methentropy(hapinfo2),ppoly(hapinfo2),length(unique(hapinfo2)), mhl(hapinfo2)),
			Fig2C=c(100*mf(hapinfo3), methentropy(hapinfo3),ppoly(hapinfo3),length(unique(hapinfo3)), mhl(hapinfo3)),
			Fig2D=c(100*mf(hapinfo4), methentropy(hapinfo4),ppoly(hapinfo4),length(unique(hapinfo4)), mhl(hapinfo4)),
			Fig2E=c(100*mf(hapinfo5), methentropy(hapinfo5),ppoly(hapinfo5),length(unique(hapinfo5)), mhl(hapinfo5)))
rownames(results)<-c("Methylation frequency", "Methylation entropy", "Epi-polymorphism", "Haplotypes", "MHL")
print(results)

# Okay-check done! -shicheng


