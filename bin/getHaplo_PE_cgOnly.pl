#!/usr/bin/perl -w
# Given SAM formated alignment file sorted by reads names, output the methylation haplotype per read
# SAM file must be sorted by read names
# Requirements: 
#   (1) CpG position table generated with genomePrep.pl (1-based)
#   (2) samtools v1.8 or above is required to be installed
# Contact: Dinh Diep (hdinhdp@gmail.com)
# Version 1.0

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


use strict;
use warnings;

my %cgTable;
my %haploTable;

my $min_coverage = 1;
my $qual_base = 33;
my $min_length = 1;
my $samtools="samtools";

my $deletion_string = "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'D'}='D';

sub main{
	printUsage() if(!$ARGV[1]);
	readCGTable($ARGV[1]);
        my $bamfile = $ARGV[0];
	my $start_time = time;
	my $TEMP_PREFIX = $ARGV[2];
	my ($last_list_pos, $last_list_calls, $last_id, $last_chr) = ("NA", "NA", "NA", "NA");
        open(INFILE, "$samtools sort -T $TEMP_PREFIX -@ 4 -n -O SAM $bamfile |") || die("Error running samtools to convert BAM file\n");
	while(my $line = <INFILE>){
		next if($line =~ m/^@/);
		chomp($line);
		my %alignment = getAlignmentInfo($line);
		getReference(\%alignment);
		my ($list_pos, $list_calls) = getCGPos(\%alignment);
		next if($list_calls eq "NA"); # skips if there is no call for this read.
		if($last_id ne $alignment{"name"} || $last_chr ne $alignment{"chr"} ){
			$haploTable{$last_list_pos}->{$last_list_calls}++ if($last_id ne "NA");
			$last_id = $alignment{"name"};
			$last_chr = $alignment{"chr"};
			$last_list_pos = $list_pos;
			$last_list_calls = $list_calls;
		}else{
			# same reads, combine list pos and list calls
			my @posLastList = split ":", $last_list_pos;
			my @lastBases = split "", $last_list_calls;
			my @posList = split ":", $list_pos;
			my @Bases = split "", $list_calls;
			my %indices;
			$last_list_pos = $last_chr;
			undef($last_list_calls);
			for(my $i = 0; $i < scalar(@Bases); $i++){
				$indices{$posList[$i+1]}=$Bases[$i];
			}
			for(my $i = 0; $i < scalar(@lastBases); $i++){
				$indices{$posLastList[$i+1]}=$lastBases[$i];
			}
			foreach my $pos (sort {$a<=>$b} keys %indices){
				$last_list_pos = $last_list_pos . ":" . $pos;
				#print $last_list_calls, "\n";
				$last_list_calls = defined($last_list_calls) ? $last_list_calls . $indices{$pos} : $indices{$pos};
				undef($indices{$pos});
			}
		
		}
	}
	$haploTable{$last_list_pos}->{$last_list_calls}++ if($last_id ne "NA");
	undef(%cgTable);
	foreach my $value(keys %haploTable){
		my @can = keys %{$haploTable{$value}};
		my @posInfo = split ":", $value;
		my ($chr, $start_pos, $end_pos) = ($posInfo[0], $posInfo[1], $posInfo[scalar(@posInfo)-1]);
		foreach my $candidate (@can){
			print "$chr:$start_pos-$end_pos", "\t", $candidate, "\t", $haploTable{$value}->{$candidate}, "\t", join(":", @posInfo[1..scalar(@posInfo)-1]), "\n";
		}
	}
	my $total_time = time - $start_time;
	warn("Total time: $total_time\n");
}

sub readCGTable{
    my $cgFile = shift;
    open(IN, "zcat $cgFile |") || die("Error opening $cgFile\n");
    while(my $line = <IN>){
        chomp($line);
        my @tmp = split /\t/, $line;
        $tmp[0] =~ s/:W//;
        $cgTable{$tmp[0].":".$tmp[1]} = 1;
    }
	close(IN);
}

sub rcSeq{	
	my $seq = shift;
	my $rcseq;
	while($seq){
		$rcseq = $rcseq . $rcTable{chop($seq)};
	}
	return $rcseq;
}

sub getReference{
	my $alignment = shift;
	my $cigar = $alignment->{"cigar"};
	my $query = $alignment->{"seq"};
	my $qual = $alignment->{"qual"};
	my $ref_seq;
	my $ref_qual;
	# Adapted from Ben Langmead (BtlBio::Alignment:Util.pm)
	# CIGAR fields are always in pairs
	my $i = 0;
	my $j = 0;
	my $nm_i = 0;
	my $nm_d = 0;
	while($i < length($cigar)){
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || die("Bad cigar string : '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || die("count not parse operation at pos $i: '$cigar'");
		die("Could not understand $op: '$cigar'") if($op !~ m/[MX=DIS]/);
		$op =~ s/[X=]/M/g;
		my ($clip_s, $clip_q);
		if($op eq "M" || $op eq "I" || $op eq "S"){
			$clip_s = substr($query, $j, $runlen);
			$clip_q = substr($qual, $j, $runlen);
			$clip_s =~ s/[ATGCatgc]/I/g if($op eq "I"); 
			$nm_i += $runlen if($op eq "I");
			$j += $runlen;
		}else{
			#deletion from reference
			$nm_d += $runlen;
			length($deletion_string) > $runlen || die("deletion is too long at $runlen: '$cigar'");
			$clip_s = substr($deletion_string, 0, $runlen);
			$clip_q = substr($deletion_string, 0, $runlen);
		}
		$i++;
		$ref_seq = $ref_seq . $clip_s if($op =~ m/[MD]/);
		$ref_qual = $ref_qual . $clip_q if($op =~ m/[MD]/);
	}
	
	$alignment->{"ref_match_seq"} = $ref_seq;
	$alignment->{"ref_match_qual"} = $ref_qual;
	$alignment->{"ref_nm_i"} = $nm_i;
	$alignment->{"ref_nm_d"} = $nm_d;
	$alignment->{"tlen"} = length($ref_seq) - $nm_i if($ref_seq);
}

sub getCGPos{
	my $alignment = shift;
	my @mismatches;
	return ("NA", "NA") if(!$alignment->{"ref_match_seq"});
	return ("NA", "NA") if(!$alignment->{"ref_match_qual"});
	my $ref_seq = uc($alignment->{"ref_match_seq"});
	my $ref_qual = $alignment->{"ref_match_qual"};
	my $nlen = length($ref_seq);
	my $chr = $alignment->{"chr"};
	my $position = $alignment->{"pos"};
	my %twomer;
	my %cgPos;
	my $actual_pos;
	$twomer{"CG"} = 1;
	if($alignment->{"orientation"} eq "C"){
		$twomer{"CA"} = 1;
		$actual_pos = $position - 1;
		if($cgTable{"$chr:$actual_pos"}){
			$b = $rcTable{substr($ref_seq,0,1)};
			$cgPos{$actual_pos} = $b if($b =~ /[CT]/);
		}
	}else{
		$twomer{"TG"} = 1;
		$actual_pos = $position + $nlen - 1;
		if($cgTable{"$chr:$actual_pos"}){
			$b = substr($ref_seq, $nlen-1, 1);
			$cgPos{$actual_pos} = $b if($b =~ /[CT]/);
		}
	}
	for(my $j = 0; $j < $nlen-1; $j++){
		my $kmer = substr($ref_seq, $j, 2);
		next if(!$twomer{$kmer});
		$actual_pos = $position + $j;
		next if(!$cgTable{"$chr:$actual_pos"});
		if($alignment->{"orientation"} eq "C"){
			$b = $rcTable{substr($kmer,1,1)};
		}else{
			$b = substr($kmer,0,1);
		}
		$cgPos{$actual_pos} = $b if($b =~ /[CT]/);
        }
	my $list_pos = $chr;
	my $list_calls;
	my @ordered_pos = sort {$a <=> $b} keys %cgPos;
	for(my $i = 0; $i < scalar(@ordered_pos); $i++){
		my $pos = $ordered_pos[$i];
		$list_pos = $list_pos . ":" . $pos;
		$list_calls = defined($list_calls) ? $list_calls . $cgPos{$pos} : $cgPos{$pos};
	}
	if(!defined($list_calls)){
		return ("NA", "NA");
	}else{
		return ($list_pos, $list_calls);
	}
}

sub getAlignmentInfo{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"name"} = (split /#/, $tmp[0])[0];
	$alignment{"orientation"} = $tmp[1] & 16 ? "C":"W";
	$alignment{"chr"} = $tmp[2];
	$alignment{"pos"} = $tmp[3];
	$alignment{"cigar"} = $tmp[5];
	$alignment{"seq"} = $tmp[9];
	$alignment{"qual"} = $tmp[10];
	#$alignment{"ndz"} = $tmp[11];
	#$alignment{"mdz"} = $tmp[12];
	$alignment{"ref_match_seq"} = "NA";
	$alignment{"ref_match_qual"} = "NA";
	$alignment{"ref_nm_i"} = 0;
	$alignment{"ref_nm_d"} = 0;
	$alignment{"tlen"} = 0;
	$alignment{"cgPos"} = ();
	return %alignment;
}

sub printUsage{
	print " Usage: \n";
	print " getHaplo_PE_cgOnly.pl [cpg position list] [sam file sorted by name] [temp_prefix]\n";
	exit 0;
}

main();
