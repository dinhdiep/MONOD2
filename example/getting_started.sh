#!/bin/bash

scripts_dir="../scripts"
bamfiles_list="cct"
metric="MHL"
target_file="small.mhbs.txt"
cpg_pos="../allcpg/hg19.fa.allcpgs.txt.gz"
outname="cct"

for bam in `cat $bamfiles_list`
do

  prefix_name=$bam
  $scripts_dir/bam2cghap.sh $cpg_pos $bam $prefix_name
  echo "$prefix_name.cgPE.hapinfo.txt" >> $outname.hap_list

done

$scripts_dir/cghap2matrix.sh $outname.hap_list $metric $outname $target_file

rm $outname.hap_list
