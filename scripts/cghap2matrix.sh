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

echo "Running $0"

files_list=$1
metric=$2
target_file=$4
outname=$3

update_files_list="$files_list"

if [ -z "$3" ] 
then
	echo "usage: $0 [list of haplotype files] <AMF|MHL|IMF> [output name prefix] [target bed file]"
	echo "          *if target bed file is not provided, assumes haplotype info do not need to be merged"
	exit 0
fi

if [ -z "$4" ]; then

	# Assume bed files have been merged
	update_files_list="$files_list"

else	
	update_files_list="$files_list.merged"
	rm $update_files_list


	# Making sure the haplotype file have been merged to regions
	for f in `cat $files_list`
	do
		bin/mergeHaploInfo_bed.pl $target_file < $f > results/HapInfo/$f.merged.hapinfo.txt
		echo "results/HapInfo/$f.merged.hapinfo.txt" >> $update_files_list
	done

fi

case "$metric" in
	AMF)
		bin/hapinfo2amf.pl $update_files_list > $outname.amf.txt
		;;
	amf)
		bin/hapinfo2amf.pl $update_files_list > $outname.amf.txt
		;;
	MHL)
		bin/hapinfo2mhl.pl $update_files_list > $outname.mhl.txt
		bin/hapinfo2umhl.pl $update_files_list > $outname.umhl.txt
		;;
	mhl)
		bin/hapinfo2mhl.pl $update_files_list > $outname.mhl.txt
		bin/hapinfo2umhl.pl $update_files_list > $outname.umhl.txt
		;;
	IMF)
		bin/hapinfo2imf.pl $update_files_list > $outname.imf.txt
		;;
	imf)
		bin/hapinfo2imf.pl $update_files_list > $outname.imf.txt

esac
