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

echo $0

name=$(basename "$0" ".sh")
src_dir=`echo $0 | sed "s/$name.sh//g" | sed "s/scripts/src/g"`
perlcode_dir="$src_dir/Perlcode"

haploInfo=$1
binFile=$2
minR2=$3
outname=$4

if [ -z "$4" ]
then
        echo "usage: $0 [haplotype file] [target bed] [minimum LD R2 cutoff] [output name prefix]"
        exit 0
fi

echo -e "#CHROM\tCHROMSTART\tCHROMEND\tBLOCKID\tNUMCPGs" > $outname.mhbs

$perlcode_dir/mergeHaploInfo_bed.pl $binFile < $haploInfo > $outname.merged.hapinfo.txt

for d in `cut -f 1 $outname.merged.hapinfo.txt | sort -u`
do	
	grep $d $outname.merged.hapinfo.txt | awk '{if(length($2) > 3) print $0}' > $outname.tmp.info
	$perlcode_dir/hapInfo2mld_block.pl $outname.tmp.info $minR2 >> $outname.mhbs
	rm $outname.tmp.info
done
