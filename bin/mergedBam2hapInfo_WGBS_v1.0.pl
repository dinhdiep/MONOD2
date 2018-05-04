#!/usr/bin/perl -w
# Bam to hapinfo files
# Requirements: 
#   (1) hg19.fa.allcpgs.txt (1-based)
#   (2) region file in BED format (chrom, start, end) 
#   (3) samtools version 1.2 or above must be installed
# Contact: Dinh Diep (hdinhdp@gmail.com)
# Version 1.0
# Created: 2015-01-14

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


die("USAGE: $0 [target bed file] [BAM file] [path to CpG positions list]\n") if(scalar(@ARGV)<3);

use strict;
my $samtools = "samtools";
my $target_list_file= $ARGV[0];
my $merged_bam_file = $ARGV[1];
my $cpg_position_file=$ARGV[2];


my $phred = 33; # in original version, 64 was used
my $target_flanking_len = 80;
my %targetTable;
my %hapInfoTable;
my %binnedtargetTable;
my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';

my ($reads_parsed, $reads_after_filtering);
my %chrSizes=("chr1","249250621","chr2","243199373","chr3","198022430","chr4","191154276",
                                "chr5","180915260","chr6","171115067","chr7","159138663",
                                "chr8","146364022","chr9","141213431","chr10","135534747",
                                "chr11","135006516","chr12","133851895","chr13","115169878",
                                "chr14","107349540","chr15","102531392","chr16","90354753",
                                "chr17","81195210","chr18","78077248","chr19","59128983",
                                "chr20","63025520","chr21","48129895", "chr22","51304566",
                                "chrX","155270560","chrY","59373566","chrM","16571");

# make a hash table for the target regions
open(INFILE, "$target_list_file")||die("Error in opening $target_list_file\n");
my $bin_size = 10000;
while(my $line = <INFILE>){
	chop($line);
	my ($chr,$target_start,$target_end,$id) = split(/[\t ]+/, $line);
	my $index = int($target_start/$bin_size);
	my $target_id = "$chr:$target_start-$target_end";
	push(@{$binnedtargetTable{$chr}->{$index}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index-1}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index+1}},$target_id);
}
close(INFILE);

# make a hash table for the CpG positions within target regions only
open(INFILE, "zcat $cpg_position_file |")||next;
while(my $line = <INFILE>){
	chop($line);
	my ($chr, $str, $pos, $context) = split(/[:\t ]+/, $line);
	$pos--; #v1.0 making sure to convert 1-base CpG positions to 0-base.
	my $index = int($pos/$bin_size);
	next if(!$binnedtargetTable{$chr}->{$index});	
	foreach my $target_id (@{$binnedtargetTable{$chr}->{$index}}){
		my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
		if($target_start <= $pos && $pos <=$target_end){
			push(@{$targetTable{$target_id}->{"CpG_positions"}},$pos);			
		}
	}
}
close(INFILE);

# for each target region, extract sequencing reads using samtools
foreach my $target_id (sort keys(%targetTable)){
	my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
	my $region_start = $target_start-$target_flanking_len;
	my $region_end = $target_end+$target_flanking_len;
	my $cmd = "$samtools view -q 10 $merged_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	next if(scalar(@lines)<1 || !$targetTable{$target_id}->{"CpG_positions"});
	$reads_parsed += scalar(@lines);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$targetTable{$target_id}->{"CpG_positions"}};
	my %readHapInfo;	

	#Going through the reads one at a time
	my %read_start_pos_table;
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $read_strand = $fields[1];		
		my $seq = $fields[9];		
		my $qual_string = $fields[10];
		my $read_start = $fields[3];
		my $CIGAR = $fields[5];
		my $read_length = length($seq);
		next if($CIGAR =~ /[ID]/);
		if($CIGAR=~ /S/){
			my ($clip_len, @others) = split(/S/, $CIGAR);
			next if(length($clip_len)>2);
			$seq=substr($seq,$clip_len,$read_length-$clip_len);		
			$qual_string=substr($qual_string,$clip_len,$read_length-$clip_len);	
			$read_length -= $clip_len;			
		}elsif($seq =~ /^CG/){
			$seq=substr($seq,2,$read_length-2);		
			$qual_string=substr($qual_string,2,$read_length-2);	
			$read_start += 2;
			$read_length -= 2;
		}
		my @read_fields = split(/[\:\#\_]+/, $fields[0]);
		my $read_id = scalar(@read_fields)>7 ? join("-", @read_fields[0..6]) : join("-", @read_fields[0..scalar(@read_fields)-1]);	
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			next if($offset<0 || $offset >= length($seq));	
			if($read_strand<16){
				my $qual_score = ord(substr($qual_string,$offset,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score);
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=substr($seq,$offset,1);
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;
			}else{
				my $qual_score = ord(substr($qual_string,$offset+1,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score);
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=$rcTable{substr($seq,$offset+1,1)};			
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;	
			}	
		}	
		$read_start_pos_table{$read_id} = $read_start if(!$read_start_pos_table{$read_id} || $read_start_pos_table{$read_id}<$read_start);
	}		
	
	#Use read start as the UMI to collapse multiple clonal reads.
	my @read_list = keys(%readHapInfo);		
	my %unique_read_base_info;
	foreach my $read_id (@read_list){		
		my $UMI = $read_start_pos_table{$read_id};
		foreach my $CpG_position (keys(%{$readHapInfo{$read_id}})){
			next if(!$readHapInfo{$read_id}->{$CpG_position}->{"base"});
			$unique_read_base_info{$UMI}->{$CpG_position}->{$readHapInfo{$read_id}->{$CpG_position}->{"base"}}+=$readHapInfo{$read_id}->{$CpG_position}->{"qual"};
		}
	}
	
	#Derive the consensus haplotype string based on multiple clonal reads.
	foreach my $UMI (keys(%unique_read_base_info)){
		my $hap_string="";
		foreach my $CpG_position(@CpG_positions){
			my $base = "N";
			my $best_qual = 0;
			if($unique_read_base_info{$UMI}->{$CpG_position}){
				foreach my $base_call (keys(%{$unique_read_base_info{$UMI}->{$CpG_position}})){
					next if($unique_read_base_info{$UMI}->{$CpG_position}->{$base_call} <= $best_qual);
					$best_qual = $unique_read_base_info{$UMI}->{$CpG_position}->{$base_call};
					$base = $base_call;
				}
			}					
			$hap_string=$hap_string.$base;
		}	
		$hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
	}

	#report haplotype strings only on the positions with valid calls.
	foreach my $hap_string(keys %{$hapInfoTable{$target_id}->{"hap_counts"}}){
		my @valid_positions;
		my $valid_hap="";
		for(my $i=0; $i<length($hap_string); $i++){
			my $allele = substr($hap_string,$i,1);
			next if($allele eq "N");
			$valid_hap = $valid_hap.$allele;
			push(@valid_positions,$CpG_positions[$i]);
		}
		next if(scalar(@valid_positions)<3);
		print $target_id,"\t$valid_hap\t", $hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@valid_positions), "\n";
	}
}
	
