#!/usr/bin/perl -w

# Hapinfo to average methylation frequency (AMF)
# Run the script to the Hapinfo directory
# Contact: Dinh Diep (hdinhdp@gmail.com)
# Version 1
# Update: 2018-03-11

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

my %coverageTable;
my %hapCountMatrix;
my @sampleList;
my $windowSize = 1000;

# making sure this region have at least 3 CpGs coverage and at least 3 reads coverage
my $min_hap_count = 3;
my $min_hap_length = 3;

sub printUsage{
  print "\nperl $0 [Regions BED] [Hap Info Table] [AMF|MHL|COV|PDR|MHL3|UMHL|UMHL3] [OUTPUT prefix] \n";
  print "Regions BED file contains the regions to generate the matrix for\n";
  print "Hap Info Table contains the table of samples to generate the matrix. First column is sampleID and second column is path to hapInfo file\n";
  print "Choose the value(s) for matrix, ie \'AMF,MHL\' will generate both AMF and MHL\n";
  print "\tAMF: average methylation frequency\n";
  print "\tMHL: methylated haplotype load\n";
  print "\tCOV: total coverage per region\n";
  print "\tPDR: percent discordant reads\n";
  print "\tMHL3: methylated haplotype load with more stringent weights\n";
  print "\tUMHL: methylated haplotype load\n";
  print "\tUMHL3: unmethylated haplotype load with more stringent weights\n";
  print "\n\n";
  exit 0;
}

sub main{
  if(!$ARGV[0]){
    printUsage();
  }
  my $region_bed = $ARGV[0];
  generateCoverageTable($region_bed);
  my $hap_files_list = $ARGV[1];

  open(INLIST, "$hap_files_list") || die("Error opening $hap_files_list\n");
  while(my $line = <INLIST>){
    chomp($line);
    my ($sample_name, $hap_file) = split "\t", $line;
    push(@sampleList, $sample_name);
    readHapTable($hap_file, $sample_name);
  }
  close(INLIST);
  
  my @valuesList = split ",", $ARGV[2];
  my $output_prefix = $ARGV[3];
  foreach my $value (@valuesList){
    generateAMF($output_prefix . ".amf.txt") if(uc($value) eq "AMF");
    generateCOV($output_prefix . ".cov.txt") if(uc($value) eq "COV");
    generatePDR($output_prefix . ".pdr.txt") if(uc($value) eq "PDR");
    generateMHL($output_prefix . ".mhl.txt", 1) if(uc($value) eq "MHL");
    generateMHL($output_prefix . ".mhl3.txt", 3) if(uc($value) eq "MHL3");
    generateUMHL($output_prefix . ".umhl.txt", 1) if(uc($value) eq "UMHL");
    generateUMHL($output_prefix . ".umhl3.txt", 3) if(uc($value) eq "UMHL3");
  }

}

sub generateAMF{

  my $output_file = shift;
  my %matrix;
  
  open(OUT, ">$output_file") || die("Error writing output $output_file\n");
  print OUT "Probe_id\t", join("\t", sort @sampleList), "\n";
  foreach my $probeID (keys(%hapCountMatrix)){
    print OUT "$probeID";
    foreach my $sample_name (sort @sampleList){
      my %k_mer_counts;
      my $mc_total=0;
      my $ct_total=0;
      foreach my $hapString (keys(%{$hapCountMatrix{$probeID}->{$sample_name}})){
        for(my $i = 0; $i < length($hapString); $i++){
          my $sub_hapString = substr($hapString,$i,1);
          next if($sub_hapString =~ /[NAG]/i);
          $mc_total+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString} if($sub_hapString eq "C");
          $ct_total+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString};
        }
      }
      if(!$ct_total || $ct_total < $min_hap_count){
         print OUT "\tNA";
      }else{
         print OUT "\t", $matrix{$probeID}->{$sample_name}=$mc_total/$ct_total;
      }
    }
    print OUT "\n";
  }
  close(OUT);
}

sub generateCOV{

  my $output_file = shift;

  open(OUT, ">$output_file") || die("Error writing output $output_file\n");
  print OUT "Probe_id\t", join("\t", sort @sampleList), "\n";
  foreach my $probeID (sort keys(%hapCountMatrix)){
    print OUT "$probeID";
    foreach my $sample_name (sort @sampleList){
      my $total_hap_counts=0;
      foreach my $hap_string (keys(%{$hapCountMatrix{$probeID}->{$sample_name}})){
        $total_hap_counts+=$hapCountMatrix{$probeID}->{$sample_name}->{$hap_string};
      }
      print OUT "\t", $total_hap_counts;
    }
    print OUT "\n";
  }
  close(OUT);
}

sub generatePDR{

  my $output_file = shift;
  
  open(OUT, ">$output_file") || die("Error writing output $output_file\n");
  print OUT "Probe_id\t", join("\t", sort @sampleList), "\n";
  foreach my $probeID (keys(%hapCountMatrix)){
    print OUT "$probeID";
    foreach my $sample_name (sort @sampleList){
      my $total_hap_counts=0;
      my $discordant_hap_counts=0;
      foreach my $hapString (keys(%{$hapCountMatrix{$probeID}->{$sample_name}})){
        $total_hap_counts+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString};
        $discordant_hap_counts+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString} if($hapString =~ m/C/i and $hapString =~ m/T/i);
      }
      if($total_hap_counts < $min_hap_count){
        print OUT "\tNA";
      }else{
        print OUT "\t", 100*$discordant_hap_counts/$total_hap_counts;
      }
    }
    print OUT "\n";
  }
  close(OUT);
}

sub generateMHL{

  my $output_file = shift;
  my $weight_exponent = shift;
  
  my @unmethylated_haps= ("T", "TT", "TTT", "TTTT", "TTTTT","TTTTTT","TTTTTTT","TTTTTTTT","TTTTTTTTT");
  my @methylated_haps  = ("C", "CC", "CCC", "CCCC", "CCCCC","CCCCCC","CCCCCCC","CCCCCCCC","CCCCCCCCC");
  
  open(OUT, ">$output_file") || die("Error writing output $output_file\n");
  print OUT "Probe_id\t", join("\t", sort @sampleList), "\n";
  foreach my $probeID (keys(%hapCountMatrix)){
    print OUT "$probeID";
    foreach my $sample_name (sort @sampleList){
      my %k_mer_counts;
      my $mc_hap_load=0;    
      my $total_hap_counts = 0;
 
      foreach my $hapString (keys(%{$hapCountMatrix{$probeID}->{$sample_name}})){
        for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
          next if($word_size>9);
          for(my $i=0; $i<=length($hapString)-$word_size; $i++){
            my $sub_hapString = substr($hapString,$i,$word_size);
            next if($sub_hapString =~ /[NAG]/i);
            $k_mer_counts{$word_size}->{$sub_hapString}+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString};          
          }
        }
        $total_hap_counts+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString};
      }

      if($total_hap_counts < $min_hap_count){
        print OUT "\tNA";
      }else{
      
        my $norm_factor=0;
        
        foreach my $word_size (keys(%k_mer_counts)){
          $k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]});
          $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]});
          my $total_count=0;
          foreach my $allele (keys(%{$k_mer_counts{$word_size}})){
            $total_count+=$k_mer_counts{$word_size}->{$allele};      
          }
          next if($total_count<1);
          my $mh_fraction = $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count;
          my $weight = $word_size**$weight_exponent;
          $mc_hap_load += $weight*$mh_fraction;
          $norm_factor+=$weight;
        }
        
        if(!$norm_factor){
          print OUT "\tNA";
        }else{
          $mc_hap_load/=$norm_factor;
          print OUT "\t$mc_hap_load";
        }
      }
    }
    print OUT "\n";
  }
  close(OUT);
}

sub generateUMHL{

  my $output_file = shift;
  my $weight_exponent = shift;
  
  ## for the sake of coding consistency, I simply switched the labels [Kun @ 2017-12-17]
  my @methylated_haps= ("T", "TT", "TTT", "TTTT", "TTTTT","TTTTTT","TTTTTTT","TTTTTTTT","TTTTTTTTT");
  my @unmethylated_haps  = ("C", "CC", "CCC", "CCCC", "CCCCC","CCCCCC","CCCCCCC","CCCCCCCC","CCCCCCCCC");
  
  open(OUT, ">$output_file") || die("Error writing output $output_file\n");
  print OUT "Probe_id\t", join("\t", sort @sampleList), "\n";
  foreach my $probeID (keys(%hapCountMatrix)){
    print OUT "$probeID";
    foreach my $sample_name (sort @sampleList){
      my %k_mer_counts;
      my $mc_hap_load=0;    
      my $total_hap_counts = 0;
 
      foreach my $hapString (keys(%{$hapCountMatrix{$probeID}->{$sample_name}})){
        for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
          next if($word_size>9);
          for(my $i=0; $i<=length($hapString)-$word_size; $i++){
            my $sub_hapString = substr($hapString,$i,$word_size);
            next if($sub_hapString =~ /[NAG]/i);
            $k_mer_counts{$word_size}->{$sub_hapString}+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString};          
          }
        }
        $total_hap_counts+=$hapCountMatrix{$probeID}->{$sample_name}->{$hapString};
      }

      if($total_hap_counts < $min_hap_count){
        print OUT "\tNA";
      }else{
      
        my $norm_factor=0;
        
        foreach my $word_size (keys(%k_mer_counts)){
          $k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]});
          $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]});
          my $total_count=0;
          foreach my $allele (keys(%{$k_mer_counts{$word_size}})){
            $total_count+=$k_mer_counts{$word_size}->{$allele};      
          }
          next if($total_count<1);
          my $mh_fraction = $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count;
          my $weight = $word_size**$weight_exponent;
          $mc_hap_load += $weight*$mh_fraction;
          $norm_factor+=$weight;
        }
        
        if(!$norm_factor){
          print OUT "\tNA";
        }else{
          $mc_hap_load/=$norm_factor;
          print OUT "\t$mc_hap_load";
        }
      }
    }
    print OUT "\n";
  }
  close(OUT);
}

sub readHapTable{
  my $hap_file = shift;
  my $sample_name = shift;
  my %haploTable;
  open(IN, $hap_file) || die("Error opening $hap_file");
  while(my $line = <IN>){
    chomp($line);
    #chr10:10000873-10001472        CCC     1       10001056:10001082:10001168
    my @tmp =  split /\t/, $line;
    my ($hap_string, $hap_count) = ($tmp[1], $tmp[2]);
    next if(scalar(@tmp) < 3);
    my ($chr, $start, $end) = split /[,:-]/, $tmp[0];
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
          $cur_hap = $cur_hap . substr($hap_string, $i, 1);
          $cur_pos = $cur_pos ? $cur_pos . ":" . $posList[$i] : $posList[$i];
        }
      }
      next if(!$cur_hap);
      $haploTable{$candidate}->{$cur_pos}->{$cur_hap} = $haploTable{$candidate}->{$cur_pos}->{$cur_hap} ? $haploTable{$candidate}->{$cur_pos}->{$cur_hap} + $hap_count : $hap_count;
      last;
    }
  }
  close(IN);
  
  foreach my $index(keys %haploTable){
  my @hapPositions = keys %{$haploTable{$index}};
  foreach my $hap_pos (@hapPositions){
    my @hapStrings = keys %{$haploTable{$index}->{$hap_pos}};
      foreach my $hap (@hapStrings){
        $hapCountMatrix{$index}->{$sample_name}->{$hap}+=$haploTable{$index}->{$hap_pos}->{$hap};
        #print $index, "\t", $hap, "\t", $haploTable{$index}->{$hapPos}->{$hap}, "\t", $hapPos, "\n";
        delete($haploTable{$index}->{$hap});
      }
      delete($haploTable{$index}->{$hap_pos});
    }
    delete($haploTable{$index});
  }
}

sub cleaningUp{
  foreach my $key (keys %coverageTable){
    delete($coverageTable{$key});
  }
  undef(%coverageTable);
}

sub generateCoverageTable{
  my $bed_file = shift;
  
  open(IN, $bed_file) || die("Error opening $bed_file");
  while(my $line = <IN>){
    chomp($line);
    my @tmp = split "[\t\ ]", $line;
    my $index = int($tmp[1]/$windowSize);
    for(my $i = $index - 5; $i <= int($tmp[2]/$windowSize) + 5 ; $i++){
      push(@{$coverageTable{$tmp[0]}->{$i}}, $tmp[0].":".$tmp[1].":".$tmp[2]);
    }
  }
  close(IN);
}


main();
