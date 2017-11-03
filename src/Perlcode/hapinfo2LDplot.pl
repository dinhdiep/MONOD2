#!/usr/bin/perl -w
use strict;
use Sort::Array qw/Sort_Table/;

my $debug=1;
my $ifplot=1;

my %hapTable; # probeID=>hapCounts=>sampleID
              #        =>CpgPositions    
			  #        =>hapCounts
			  # 	   =>totalHap

my %LD_matrix;
my %COV_matrix;

my ($region_chr, $region_start, $region_end) = split /[:-]/, $ARGV[1];
			  
sub main(){
  my $hapInfo_file = $ARGV[0];
  load_hapInfo_files($hapInfo_file); # here filtered all regions outside of target 
  my @all_probe_IDs = sort keys(%hapTable);

  my $outname = "$region_chr,$region_start,$region_end";
  open(RSQ_OUT, ">$outname.rsq") || die("error writing to file $outname.rsq\n");
  open(COV_OUT, ">$outname.cov") || die("error writing to file $outname.cov\n");

  foreach my $probeID (@all_probe_IDs){
    my ($chr, $region_start, $region_end) = split(/[\-\:]/, $probeID);
    my @unique_cpg_positions = @{$hapTable{$probeID}->{"cpgPositions"}};
    for(my $i = 0; $i < scalar(@unique_cpg_positions) - 1; $i++){
      for(my $j = $i+1; $j < scalar(@unique_cpg_positions); $j++){
        my $rsq = lookupLD($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $i, $j);
        $LD_matrix{$unique_cpg_positions[$i]}->{$unique_cpg_positions[$j]} = $rsq;
        $LD_matrix{$unique_cpg_positions[$j]}->{$unique_cpg_positions[$i]} = $rsq;
      }
    }
  }
  my @sorted_list = sort {$a<=>$b} keys %COV_matrix;
  print RSQ_OUT "\t", join("\t", @sorted_list), "\n";
  for(my $i = 0; $i < scalar(@sorted_list); $i++){
    print RSQ_OUT $sorted_list[$i];
    print COV_OUT $sorted_list[$i], ",", $COV_matrix{$sorted_list[$i]}, "\n";
    for(my $j = 0; $j < scalar(@sorted_list); $j++){
      if(!$LD_matrix{$sorted_list[$i]}->{$sorted_list[$j]}){
        print RSQ_OUT "\tNA";
        next;
      }
      print RSQ_OUT "\t", $LD_matrix{$sorted_list[$i]}->{$sorted_list[$j]};
    }
    print RSQ_OUT "\n";
  }
	
  close(RSQ_OUT);
  close(COV_OUT);
  if($ifplot){
    plotRSQ($outname, "$outname.rsq", "$outname.cov");
  }
}			


sub plotRSQ {
   
   my $outname = shift;
   my $rsq_file = shift;
   my $cov_file = shift;
   
   open(RSCRIPT, ">drawLDheatmap.R") || die("Error writing draw heatmap R script\n");
   print RSCRIPT "\n
     args = commandArgs(trailingOnly=T)\n
     INFILE = args[1]\n
     COVERAGES = args[2]\n
     OUTFILE = args[3]\n
     mycols = c(1,   1,   1,   2,   3,   4,   5,   5,   5,  5)\n
     library(LDheatmap)\n
     mydata <- read.table(INFILE, header=T, row=1, sep=\"\\t\")\n
     mydata.transformed <- as.matrix(mydata)\n
     mydata.transformed[which(mydata.transformed < 0.3)] = 0\n
     mydata.transformed[which(mydata.transformed > 0.2 & mydata.transformed < 0.4)] = 0.2\n
     mydata.transformed[which(mydata.transformed >= 0.4 & mydata.transformed < 0.5)] = 0.4\n
     mydata.transformed[which(mydata.transformed >= 0.5 & mydata.transformed < 0.6)] = 0.6\n
     mydata.transformed[which(mydata.transformed >= 0.6 & mydata.transformed < 0.8)] = 0.8\n
     mydata.transformed[which(mydata.transformed >= 0.8 & mydata.transformed <= 1)] = 1\n
     colnames(mydata.transformed) <- rownames(mydata)\n
     rownames(mydata.transformed) <- rownames(mydata)\n
     mycoverages <- read.table(COVERAGES, header=F, row=1, sep=\",\",stringsAsFactors=F)\n
     print(mycoverages)\n
     pdf(paste(OUTFILE,\".pdf\", sep=\"\"), height=3, width=6)\n
     ld.map<-LDheatmap(mydata.transformed, title=NULL, SNP.name=NULL, flip=TRUE, color=heat.colors(5), genetic.distances=as.numeric(rownames(mydata)))\n
     LDheatmap.addScatterplot(ld.map, as.matrix(mycoverages), type=\"lines\")\n
     dev.off()\n";
  close(RSCRIPT);
  
  my $cmd = "Rscript drawLDheatmap.R $rsq_file $cov_file $outname";
  system($cmd);

  # cleaning up
  unlink("drawLDheatmap.R");
  unlink($rsq_file);
  unlink($cov_file);
  
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
		my ($chr, $start, $end) = split /[:,-]/, $probeID;
		my $hapString = $fields[1];
		next if(length($hapString)<3 || $chr ne $region_chr);
		my $hapCount = $fields[2];
		my @cpgPositions = split(/[,:]/, $fields[3]);
		push(@{$hapTable{$probeID}->{"hapInfo"}->{"Sample"}}, "$hapCount:$hapString:".$fields[3]);
		foreach my $pos(@cpgPositions){
			next if($pos < $region_start || $pos > $region_end);
			$COV_matrix{$pos} += $fields[2];
			$hapTable{$probeID}->{"cpgPositionTable"}->{$pos}=1;
		}
	}
	close(INFILE);
	
	my @all_probe_IDs = sort keys(%hapTable);
	foreach my $probeID (@all_probe_IDs){
		if(	!$hapTable{$probeID}->{"hapInfo"}->{"Sample"} ){
			delete($hapTable{$probeID});
			next;
		}

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
				next if(! $COV_matrix{$pos} );
				$allele_list[$pos2index{$pos}]=substr($hapString,$i,1);
			}
			my $full_hap_string = join("", @allele_list);
			$hapTable{$probeID}->{"hapCounts"}->{"Sample"}->{$full_hap_string}+=$count;
			$hapTable{$probeID}->{"totalHap"}->{"Sample"}+=$count;
		}		
		
		#print "$probeID\t", join(",", @unique_pos), "\n";
		delete($hapTable{$probeID}->{"hapInfo"});
		delete($hapTable{$probeID}->{"cpgPositionTable"});
		undef(%pos2index);		
	}	
}


sub lookupLD(){
	my ($h_hap_count_table, $locusA, $locusB) = @_;	
	
	my %two_locus_hapTable;
		
	$two_locus_hapTable{"CC"}=0;
	$two_locus_hapTable{"TT"}=0;
	$two_locus_hapTable{"CT"}=0;
	$two_locus_hapTable{"TC"}=0;

	foreach my $hapString (keys(%{$h_hap_count_table})){
		my $two_locus_hap = substr($hapString,$locusA,1).substr($hapString,$locusB,1);
		next if($two_locus_hap =~ /[AGN]/);
		$two_locus_hapTable{$two_locus_hap} += ${$h_hap_count_table}{$hapString};
	}
	
	my ($abs_Dprime,$r2, $abs_d, $abs_Q) = haplotype2LD(\%two_locus_hapTable);
	
	$r2=-0.1 if($r2 =~ /NA/);
	return $r2;
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
