#!/usr/bin/perl -w
# Hapinfo to methylation haplotype load (MHL)
# Run the script to the Hapinfo directory
# Contact: Kun Zhang
# Version 1.3
# Update: 2016-02-29

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
use Cwd;
die &USAGE if @ARGV <1;
my %mch_load_matrix;
my %probe_HMH_samples;
my %hap_count_matrix;
my $hapinfList=shift @ARGV;
open FF,$hapinfList;
chomp(my @hapInfo_files=<FF>);
close FF;

my @sample_list;
foreach my $hapInfo_file(sort @hapInfo_files){
        my @line=split /\//,$hapInfo_file;
	my $sample_name = "fileID_".$line[$#line];
	$sample_name =~ s/.hapInfo.txt//;
	push(@sample_list, $sample_name);
	open(INFILE, "$hapInfo_file") || die("Error in opening $hapInfo_file!");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<1);		
		my @cgPos = split /[:,]/, $fields[3];
		for(my $i = 0; $i < length($hapString); $i++){
			my $cur_pos = $cgPos[$i];
			my $sub_hapString = substr($hapString,$i,1);
			next if($sub_hapString =~ /[NAG]/i);
			$mch_load_matrix{$probeID.":".$cur_pos}->{$sample_name}->{"mc_total"}+=$fields[2] if($sub_hapString eq "C");
			$mch_load_matrix{$probeID.":".$cur_pos}->{$sample_name}->{"ct_total"}+=$fields[2];
		}
	}
	close(INFILE);
}


print "Probe_id\t", join("\t", sort @sample_list), "\n";
foreach my $probeID (sort keys(%mch_load_matrix)){
	print "$probeID";
	foreach my $sample_name(sort @sample_list){
		if(!$mch_load_matrix{$probeID}->{$sample_name}->{"ct_total"}){
			print "\tNA";
		}else{
			$mch_load_matrix{$probeID}->{$sample_name}->{"mc_total"} = 0 if(!$mch_load_matrix{$probeID}->{$sample_name}->{"mc_total"});
			my $mf = $mch_load_matrix{$probeID}->{$sample_name}->{"mc_total"}/$mch_load_matrix{$probeID}->{$sample_name}->{"ct_total"};
			print "\t", sprintf("%4.3f", $mf);
		}
	}
	print "\n";
}

sub USAGE{
	print "\nperl $0 Hapinfo_File_list > Ouput.txt\n";
	print "Just use: ls *hapInfo.txt > Hapinfo_File_list to Get Hapinfo_File_list\n";
}
