#!/usr/bin/perl -w

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

my %haploTable;
my %index2PosList;
my %probeID;
my %coverageTable;
my $windowSize = 1000;
my $minLength = 1;

# remove haplotypes with improperly merged PE reads

open(IN, $ARGV[0]) || die("Error opening $ARGV[0]");
while(my $line = <IN>){
     chomp($line);
     my @tmp = split "\t", $line;
     my $index = int($tmp[1]/$windowSize);
     for(my $i = $index - 5; $i <= int($tmp[2]/$windowSize) + 5 ; $i++){
	push(@{$coverageTable{$tmp[0]}->{$i}}, $tmp[0].":".$tmp[1].":".$tmp[2]);
     }
}
close(IN);

while(my $line = <STDIN>){
	chomp($line);
	#chr10:10000873-10001472        CCC     1       10001056:10001082:10001168
	my @tmp =  split /\t/, $line;
	my ($hapString, $hapCount) = ($tmp[1], $tmp[2]);
	next if(length($hapString) < $minLength);
	my ($chr, $start, $end) = split /[:-]/, $tmp[0];
	$tmp[3] =~ s/,/:/g;
	my @posList = split ":", $tmp[3];
	$start = $posList[0];
	$end = $posList[scalar(@posList)-1];
	next if($end - $start > 400); # discard paired reads that are too far apart
	my $bin = int($posList[0]/$windowSize);
	next if(!$coverageTable{$chr}->{$bin});
	my @can = @{$coverageTable{$chr}->{$bin}};
	foreach my $candidate(@can){
		my ($q_chr, $q_start, $q_end) = split /[:-]/, $candidate;
		next if($q_chr ne $chr || $start > $q_end || $end < $q_start);
		# put the reads in this bin
		my $cur_hap="";
		my $cur_pos="";
		for(my $i = 0; $i < scalar(@posList); $i++){
			if($posList[$i] >= $q_start and $posList[$i] <= $q_end){
				$cur_hap = $cur_hap . substr($hapString, $i, 1);
				$cur_pos = $cur_pos ? $cur_pos . ":" . $posList[$i] : $posList[$i];
			}
		}
		next if(!$cur_hap);
		$haploTable{$candidate}->{$cur_pos}->{$cur_hap} = $haploTable{$candidate}->{$cur_pos}->{$cur_hap} ? $haploTable{$candidate}->{$cur_pos}->{$cur_hap} + $hapCount : $hapCount;
		last;
	}
}

foreach my $index(keys %haploTable){
	my @hapPositions = keys %{$haploTable{$index}};
	foreach my $hapPos (@hapPositions){
		my @hapStrings = keys %{$haploTable{$index}->{$hapPos}};
		foreach my $hap (@hapStrings){
			print $index, "\t", $hap, "\t", $haploTable{$index}->{$hapPos}->{$hap}, "\t", $hapPos, "\n";
			delete($haploTable{$index}->{$hap});
		}
		delete($haploTable{$index}->{$hapPos});
	}
	delete($haploTable{$index});
}

foreach my $key (keys %coverageTable){
	delete($coverageTable{$key});
}
undef(%haploTable);
undef(%coverageTable);
