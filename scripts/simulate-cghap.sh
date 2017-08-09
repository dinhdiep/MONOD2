#!/bin/bash

# This software is Copyright © 2017 The Regents of the University of California. All Rights Reserved.
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# 858) 534-5815
# invent@ucsd.edu
# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


# Make simulations of fixed sizes at ratios 1:5, 1:10, 1:20, 1:100, 0:1

#################################################################################
################# Modify the following paths and variables ######################
#################################################################################

numsimulation=2
n1=1000

cctfastq="BAMfiles/CCT.readIDs.txt"
lctfastq="BAMfiles/LCT.readIDs.txt"
ncpfastq="BAMfiles/NCP.readIDs.txt"

cctbam="BAMfiles/CCT.bam"
lctbam="BAMfiles/LCT.bam"
ncpbam="BAMfiles/NCP.bam"

allcpgs="allcpg/hg19.fa.allcpgs.txt.gz"


################ DO NOT MODIFY BELOW THIS LINE ##################################

name=$(basename "$0" ".sh")
src_dir=`echo $0 | sed "s/$name.sh//g" | sed "s/scripts/src/g"`
scripts_dir=`echo $0 | sed "s/$name.sh//g"`
perlcode_dir="$src_dir/Perlcode"

mkdir -p "Simulated"
mkdir -p "LCT_simulation"
mkdir -p "CCT_simulation"
mkdir -p "NCP_simulation"

echo -e "minFraction\tforeground\tbackground" > cct.hapinfo.list_$numsimulation
echo -e "minFraction\tforeground\tbackground" > lct.hapinfo.list_$numsimulation
echo "" > NCP_simulation/empty.hapinfo.txt

zcat $allcpgs > Simulated/allcpgs

cpgs="Simulated/allcpgs"

for i in `seq 1 1 $numsimulation`
do

  outprefix="NCP.$RANDOM"
  $perlcode_dir/sample_reads_from_fqnames.pl $ncpfastq $ncpbam $n1 > Simulated/${outprefix}.sam
  samtools view -Sb -t BAMfiles/hg19.fa.fai Simulated/${outprefix}.sam > Simulated/${outprefix}.bam
  $scripts_dir/bam2cghap.sh $cpgs Simulated/${outprefix}.bam NCP_simulation/${outprefix}

  echo -e "0.00\tNCP_simulation/${outprefix}.bam.cgPE.hapinfo.txt\tNCP_simulation/empty.hapinfo.txt" >> cct.hapinfo.list_$numsimulation
  echo -e "0.00\tNCP_simulation/${outprefix}.bam.cgPE.hapinfo.txt\tNCP_simulation/empty.hapinfo.txt" >> lct.hapinfo.list_$numsimulation

  for x in 0.01 0.05 0.10 0.20 
  do
    n2=`echo "$n1 - $x * $n1" | bc | awk '{printf("%d", $1)}'`
    n3=`expr $n1 - $n2`
  
    outprefix1="NCP.$RANDOM"
    outprefix2="LCT.$RANDOM"
    outprefix3="CCT.$RANDOM"

    echo -e "${x}\tNCP_simulation/${outprefix1}.cgPE.hapinfo.txt\tCCT_simulation/${outprefix3}.cgPE.hapinfo.txt" >> cct.hapinfo.list_$numsimulation
    echo -e "${x}\tNCP_simulation/${outprefix1}.cgPE.hapinfo.txt\tLCT_simulation/${outprefix2}.cgPE.hapinfo.txt" >> lct.hapinfo.list_$numsimulation

    $perlcode_dir/sample_reads_from_fqnames.pl $ncpfastq $ncpbam $n2 > Simulated/${outprefix1}.sam
    samtools view -Sb -t BAMfiles/hg19.fa.fai Simulated/${outprefix1}.sam > Simulated/${outprefix1}.bam
    $scripts_dir/bam2cghap.sh $cpgs Simulated/${outprefix1}.bam NCP_simulation/${outprefix1}

    $perlcode_dir/sample_reads_from_fqnames.pl $cctfastq $cctbam $n3 > Simulated/${outprefix2}.sam
    samtools view -Sb -t BAMfiles/hg19.fa.fai Simulated/${outprefix2}.sam > Simulated/${outprefix2}.bam
    $scripts_dir/bam2cghap.sh $cpgs Simulated/${outprefix2}.bam CCT_simulation/${outprefix2}

    $perlcode_dir/sample_reads_from_fqnames.pl $lctfastq $lctbam $n3 > Simulated/${outprefix3}.sam
    samtools view -Sb -t BAMfiles/hg19.fa.fai Simulated/${outprefix3}.sam > Simulated/${outprefix3}.bam
    $scripts_dir/bam2cghap.sh $cpgs Simulated/${outprefix3}.bam LCT_simulation/${outprefix3}

  done

done

rm -r Simulated
