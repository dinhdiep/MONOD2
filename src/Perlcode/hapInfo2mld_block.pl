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
use Sort::Array qw/Sort_Table/;

my %hapTable; # probeID=>hapCounts=>sampleID
              #        =>CpgPositions    
			  #        =>hapCounts
			  # 	   =>totalHap

my %LD_matrix;			  
sub main(){
	my $hapInfo_file = $ARGV[0];
	my $min_r2 = $ARGV[1];
	$min_r2 = 0.1 if (!$min_r2);
	load_hapInfo_files($hapInfo_file);	
	my @all_probe_IDs = sort keys(%hapTable);
	foreach my $probeID (@all_probe_IDs){
		undef(%LD_matrix);
		get_methylation_LD_blocks_greedy($probeID,$min_r2);

		my @unique_cpg_positions = @{$hapTable{$probeID}->{"cpgPositions"}};
		foreach my $block_id (sort keys(%{$hapTable{$probeID}->{'MLD_blocks'}})){
			my $n_cpg_in_block = $hapTable{$probeID}->{'MLD_blocks'}->{$block_id}->{'end'} - $hapTable{$probeID}->{'MLD_blocks'}->{$block_id}->{'start'}+1; # make sure the size is +1
			next if($n_cpg_in_block<3);
			my ($chr, $region_start, $region_end) = split(/[\-\:]/, $probeID);
			print $chr,"\t", 
				$unique_cpg_positions[$hapTable{$probeID}->{'MLD_blocks'}->{$block_id}->{'start'}], "\t",
				$unique_cpg_positions[$hapTable{$probeID}->{'MLD_blocks'}->{$block_id}->{'end'}]+1, "\t", # make sure the last position is plus 1
				"$probeID,$block_id\t$n_cpg_in_block\n";
		}
	}
}			

  
sub load_hapInfo_files(){
	my $hapInfo_file = shift;
	open(INFILE, "$hapInfo_file")||die("Error in opening file $hapInfo_file\n");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		$fields[3] =~ s/:/,/g;
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<3);
		my $hapCount = $fields[2];
		my @cpgPositions = split(/[,:]/, $fields[3]);
		push(@{$hapTable{$probeID}->{"hapInfo"}->{"Sample"}}, "$hapCount:$hapString:".$fields[3]);
		foreach my $pos(@cpgPositions){
			$hapTable{$probeID}->{"cpgPositionTable"}->{$pos}=1;
		}
	}
	close(INFILE);
	#print "$hap_number haplotypes loaded for $hapInfo_file\n";
	
	my @all_probe_IDs = sort keys(%hapTable);
	foreach my $probeID (@all_probe_IDs){
		if(	!$hapTable{$probeID}->{"hapInfo"}->{"Sample"} ){
			delete($hapTable{$probeID});
			next;
		}
		my ($chr, $target_start, $target_end) = split(/[:\-]/, $probeID);
		my %cpg_pos_table;
		my %pos2index;
		
		
		my @unique_pos = sort {$a <=> $b} keys(%{$hapTable{$probeID}->{"cpgPositionTable"}});
		@{$hapTable{$probeID}->{"cpgPositions"}} = @unique_pos;
		my $full_hap_length = scalar(@unique_pos);
		for(my $i=0; $i<$full_hap_length; $i++){
			$pos2index{$unique_pos[$i]}=$i;
		}
		
		foreach my $this_hap_info (@{$hapTable{$probeID}->{"hapInfo"}->{"Sample"}}){
			my ($count, $hapString, $pos_string) = split(/:/, $this_hap_info);
			my @this_hap_positions = split(/,/, $pos_string);
			my @allele_list;
			for(my $i=0; $i<$full_hap_length; $i++){
				push(@allele_list,"N");
			}
			for(my $i=0; $i<scalar(@this_hap_positions); $i++){
				my $pos=$this_hap_positions[$i];
				$allele_list[$pos2index{$pos}]=substr($hapString,$i,1);
				$cpg_pos_table{$pos}+=$count if(substr($hapString,$i,1) =~ /[CT]/);
			}
			my $full_hap_string = join("", @allele_list);
			$hapTable{$probeID}->{"hapCounts"}->{"Sample"}->{$full_hap_string}+=$count;
			$hapTable{$probeID}->{"totalHap"}->{"Sample"}+=$count;
		}		
		
		#print "$probeID\t", join(",", @unique_pos), "\n";
		delete($hapTable{$probeID}->{"hapInfo"});
		delete($hapTable{$probeID}->{"cpgPositionTable"});
		undef(%cpg_pos_table);
		undef(%pos2index);		
	}	
}

sub get_methylation_LD_blocks_greedy(){
	my $probeID = shift;
	my $threshold = shift;
	my $n_loci = scalar(@{$hapTable{$probeID}->{"cpgPositions"}});
	
	my @blocks;
	my %blockSetTable;
	my $start = 0;
	my $i = 0;
	while($i <= $n_loci){
		while(($i < $n_loci) && (&lookupLD($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $i, $i+1) > $threshold)){
			$i++;
		}
		#print "Search between $start-$i\n";
		push(@blocks, getAllBlockInRegion($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $start, $i, $threshold));
		$i++;
		$start = $i;
	}
	my @sortedBlocks = Sort_Table(
			cols		=> '2',
	        field		=> '1',
			sorting  	=> 'ascending',
			structure	=> 'csv',
			separator	=> '\:',
			data		=> \@blocks,
	);
	
	for(my $i=0; $i<scalar(@blocks); $i++){
		my $id = sprintf("B%03d", $i+1);
		my @words = split(/:/, $sortedBlocks[$i]);
		$hapTable{$probeID}->{'MLD_blocks'}->{$id}->{'start'} = $words[0];
		$hapTable{$probeID}->{'MLD_blocks'}->{$id}->{'end'} = $words[1];
	}
}


sub getAllBlockInRegion(){
	my $h_hap_count_table = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	#print "Cur $start,$end,$threshold\n";
	my @blocks;
	if($start == $end){
		push(@blocks,$start . ":" . $end);
		return @blocks;
	}
	my ($block_start, $block_end) = findMaxBlockInRegion($h_hap_count_table, $start, $end, $threshold);
	push(@blocks, $block_start. ":" . $block_end);
	if($block_start > $start){
		my @sub_blocks = &getAllBlockInRegion($h_hap_count_table, $start, $block_start-1, $threshold);
		push(@blocks, @sub_blocks);
	}
	if($block_end < $end){
		my @sub_blocks = &getAllBlockInRegion($h_hap_count_table, $block_end+1, $end, $threshold);
		push(@blocks, @sub_blocks);
	}
	return @blocks;
}

sub findMaxBlockInRegion(){
	my $h_hap_count_table = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my $max_block_start=$start;
	my $max_block_end=$start;
	return ($start, $end) if($start == $end);
	for(my $size = $end-$start+1; $size >1; $size--){
		my $good_block=0;
		for(my $i= $start; $i<= $end-$size+1; $i++){
			$good_block = 1;
			for(my $j= $i; $j<$i+$size; $j++){
				for(my $k = $j+1; $k<$i+$size; $k++){
					if((abs($j-$k)<5) # don't check CpG sites that are too far apart
						&& lookupLD($h_hap_count_table, $j, $k) < $threshold) {
						$good_block = 0;
						last;
					}
				}
			}
			if($good_block){
				$max_block_start = $i;
				$max_block_end = $i + $size -1;
				last;
			}
		}
		last if($good_block);
	}
	return ($max_block_start,  $max_block_end);
}

sub lookupLD(){
	my ($h_hap_count_table, $locusA, $locusB) = @_;	
	if (!$LD_matrix{$locusA}->{$locusB}){	
		my %two_locus_hapTable;
		#$two_locus_hapTable{"00"}=0;
		#$two_locus_hapTable{"01"}=0;
		#$two_locus_hapTable{"11"}=0;
		#$two_locus_hapTable{"10"}=0;
		
		$two_locus_hapTable{"CC"}=0;
		$two_locus_hapTable{"TT"}=0;
		$two_locus_hapTable{"CT"}=0;
		$two_locus_hapTable{"TC"}=0;

		foreach my $hapString (keys(%{$h_hap_count_table})){
			my $two_locus_hap = substr($hapString,$locusA,1).substr($hapString,$locusB,1);
			next if($two_locus_hap =~ /[AGN]/);
			#$two_locus_hap =~ s/C/1/g;
			#$two_locus_hap =~ s/T/0/g;
			$two_locus_hapTable{$two_locus_hap} += ${$h_hap_count_table}{$hapString};
		}
	
		my ($abs_Dprime,$r2, $abs_d, $abs_Q) = haplotype2LD(\%two_locus_hapTable);
	
		#print "#$locusA,$locusB:", "\tCC=",$two_locus_hapTable{"CC"},  "\tCT=",$two_locus_hapTable{"CT"}, "\tTC=",$two_locus_hapTable{"TC"},"\tTT=",$two_locus_hapTable{"TT"},"\tr2=$r2\n";
		
		$r2=-0.1 if($r2 =~ /NA/);
		$r2+=0.0001;
		$LD_matrix{$locusA}->{$locusB} = $r2;
		$LD_matrix{$locusB}->{$locusA} = $r2;
	}
	#print "$locusA,$locusB\t r2=",$LD_matrix{$locusA}->{$locusB},"\n";
	return $LD_matrix{$locusA}->{$locusB}-0.0001;
}


sub haplotype2LD(){
	my $h_two_locus_hapTable = shift;
	my %two_locus_hapTable = %{$h_two_locus_hapTable};
	my @AlleleTable;
	my ($abs_Dprime,$r2, $abs_d, $abs_Q, $mafA, $mafB);
	foreach my $hap (keys(%two_locus_hapTable)){
		if(!$AlleleTable[0][0]){
				$AlleleTable[0][0] = substr($hap,0,1);
		}elsif($AlleleTable[0][0] ne substr($hap,0,1)){
				$AlleleTable[0][1] = substr($hap,0,1);
		}
		if(!$AlleleTable[1][0]){
				$AlleleTable[1][0] = substr($hap,1,1);
		}elsif($AlleleTable[1][0] ne substr($hap,1,1)){
				$AlleleTable[1][1] = substr($hap,1,1);
		}
	}

	foreach my $alleleA (@{$AlleleTable[0]}){
		foreach my $alleleB (@{$AlleleTable[1]}){
			if(!$two_locus_hapTable{$alleleA.$alleleB}){
				$two_locus_hapTable{$alleleA.$alleleB} = 0.0;
			}
		}
	}

	my @AlleleFreq;
	$AlleleFreq[0][0] = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} + $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]};
	$AlleleFreq[0][1] = $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]};
	$AlleleFreq[1][0] = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]};
	$AlleleFreq[1][1] = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]} + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]};
	return ("NA","NA","NA","NA", "NA", "NA") if($AlleleFreq[0][0] + $AlleleFreq[0][1] ==0 || $AlleleFreq[1][0] + $AlleleFreq[1][1] == 0);

	$mafA = $AlleleFreq[0][0] < $AlleleFreq[0][1] ? $AlleleFreq[0][0] : $AlleleFreq[0][1];
	$mafA /= $AlleleFreq[0][0] + $AlleleFreq[0][1];
	$mafB = $AlleleFreq[1][0] < $AlleleFreq[1][1] ? $AlleleFreq[1][0] : $AlleleFreq[1][1];
	$mafB /= $AlleleFreq[1][0] + $AlleleFreq[1][1];
	return ("NA","NA","NA","NA", $mafA, $mafB) if($mafA ==0 || $mafB==0);
	
	my $D = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]}
		- $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]};
	my $Dmax = $D > 0 ? ($AlleleFreq[0][0]*$AlleleFreq[1][1] < $AlleleFreq[0][1]*$AlleleFreq[1][0] ? $AlleleFreq[0][0]*$AlleleFreq[1][1]:$AlleleFreq[0][1]*$AlleleFreq[1][0]) :
				($AlleleFreq[0][0]*$AlleleFreq[1][0] < $AlleleFreq[0][1]*$AlleleFreq[1][1] ? $AlleleFreq[0][0]*$AlleleFreq[1][0] : $AlleleFreq[0][1]*$AlleleFreq[1][1]);

	if($D == 0.0 ){
		$abs_Dprime = ($Dmax == 0.0)? 1.0 :0.0;
		$r2 = ($AlleleFreq[0][0]*$AlleleFreq[0][1]*$AlleleFreq[1][0]*$AlleleFreq[1][1]) == 0.0 ?
				1.0 : 0.0;
		$abs_d = $AlleleFreq[1][0]*$AlleleFreq[1][1] == 0.0 ? 1.0 : 0.0;
		$abs_Q = ($two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]}
				+ $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]} == 0.0) ?
				1.0 : 0.0;
	}else{
		$abs_Dprime = abs($D/$Dmax);
		$r2 = $D*$D/($AlleleFreq[0][0]*$AlleleFreq[0][1]*$AlleleFreq[1][0]*$AlleleFreq[1][1]);
		$abs_d = abs($D/($AlleleFreq[1][0]*$AlleleFreq[1][1]));
		$abs_Q = abs($D/($two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]}
		     + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]}));
	}
	return ($abs_Dprime,$r2, $abs_d, $abs_Q, $mafA, $mafB);
}

main();
