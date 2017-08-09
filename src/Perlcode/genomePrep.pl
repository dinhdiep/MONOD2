#!/usr/bin/perl -w
# Usage 1: genomePrep.pl genome.fa[.gz] 
# convert: no = don't convert genome
# context: CG = CG only
# Requirements: 
#   (1) genome fasta file
#   (2) gunzip required for gzip fasta files
# Contact: Dinh Diep (hdinhdp@gmail.com)
# Version 1.0

use strict;

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';

my %variant;
my $context = 1;
my $convert = 0;
my $all = 0;

my $genome_path = $ARGV[0];
printUsage("No genome provided.\n") if(!$ARGV[0]);
my @tmp = split /\//, $genome_path;
my $genome = pop(@tmp);
undef @tmp;


for(my $i = 1; $i<scalar(@ARGV); $i++){
	my $val = $ARGV[$i];
	if($val =~ m/convert/){ $val=~s/convert=//g; $convert = 0 if($val =~ m/no/i); }
	if($val =~ m/context/){ $context = 1; $val=~s/context=//g; $all = 1 if($val =~ m/all/i); $context = 0 if($val =~ m/none/i);}
}

my $bisCT = $genome . ".bis.CT";
my $bisGA = $genome . ".bis.GA";

sub main{
	my $start = time;
	my $variant_file = $ARGV[1];
	read_variants($variant_file) if($variant_file);
	if($convert){
	        open(FWD_OUT, ">$bisCT") || die("Error writing C->T file");
	        close(FWD_OUT);
	        open(REV_OUT, ">$bisGA") || die("Error writing G->A file");
        	close(REV_OUT);
	}
	process_genome($genome_path);
	my $time_taken = time - $start;
}

sub read_variants{
	my $variant_file = shift;
	if(open(IN, "$variant_file")){

	}else{
		print "No variant file provided or variant file is not readable\n";
		print "Continuing without variant file\n";
		return 0;
	}
	$genome = "var." . $genome;
	while(my $line = <IN>){
		chop($line);
		my @f = split /\t/, $line;
		next if($f[10] ne "ref" && $f[16] ne "ref");
		$variant{$f[2]}->{$f[4]} = $f[11] if($f[10] eq "snp");
		$variant{$f[2]}->{$f[4]} = $f[17] if($f[16] eq "snp");
		#print $f[3], "\t", $f[5], "\t",  $f[10], "\t", $f[16], "\t", $f[11], "\t", $f[17], "\n";
	}
	close(IN);
}

sub process_genome{
	my $genome = shift;
	if($genome =~ /\.gz$/) {
		open(IN, "gunzip -c $genome |") || printUsage("Can't open pipe to $genome\n");
	}else{
		open(IN, $genome) || printUsage("Can't open $genome\n");
	}
	my $line = <IN>;
	chop($line);
	while($line =~ m/>/){
		$line =~ s/>//g;	
		my $name = $line;
		#die if($name ne "chr10"); ###remember to remove!!
		my $fwd_seq = "NA";
		my @rcArray;
		my $pos = 0;
		my ($first, $second, $third) = ('N', 'N', 'N');
		if($context){
			open(CH_OUT, ">$genome.$name.cpositions.txt") || die("Error writing c positions file");
		}
		if($convert){
			open(CT_OUT, ">>$bisCT") || die("Error writing C->T file");
			print CT_OUT ">", $name, "_Watson\n";
			open(GA_OUT, ">>$bisGA") || die("Error writing G->A file");
			print GA_OUT ">", $name, "_Crick\n";
		}
		while(my $seq = <IN>){	
			chop($seq);
			$seq =~ tr/ //;
			next if(!$seq);
			if($seq =~ m/>/){
				$line = $seq;
				last;
			}
			$seq = uc($seq);
			my @base = split("", $seq);
			for(my $i = 0; $i < scalar(@base); $i++){
				$pos++;
				my $val = $base[$i];
				if($variant{$name}->{$pos}){
					$val = $variant{$name}->{$pos};
					$base[$i] = $val;
				}
				if($convert){
					if($fwd_seq eq "NA"){
						$fwd_seq = $val;
					}else{
						$fwd_seq = $fwd_seq . $val;
					}
				}
				if($pos > 3 && $context){
					my $c_pos = $pos - 3;
					print CH_OUT $name, ":W\t", $c_pos, "\tCG\n" if($first eq 'C' and $second eq 'G');
					if($all){
						if($first eq 'C' and $second ne 'G'){
							print CH_OUT $name, ":W\t", $c_pos, "\tCHG\n", $name, ":C\t", $c_pos+2, "\tCHG\n" 
								if($third eq 'G' and $second ne 'C');
							print CH_OUT $name, ":W\t", $c_pos, "\tCHG\n"
								if($third eq 'G' and $second eq 'C');
							print CH_OUT $name, ":W\t", $c_pos, "\tCHH\n" 
								if($third ne 'G');
						}
						print CH_OUT $name, ":C\t", $c_pos+2, "\tCHH\n" 
							if($first ne 'C' and $second ne 'C' and $third eq 'G');
						print CH_OUT $name, ":C\t", $c_pos+2, "\tCHG\n"
							if($first eq 'C' and $second eq 'G' and $third eq 'G');
					}
				}
				($first, $second, $third) = ($second, $third, $val);
				if($convert && $pos%150 == 0){
					my ($ct_seq, $ga_seq) = ($fwd_seq, $fwd_seq);
					$ct_seq =~ tr/C/T/;
					$ga_seq =~ tr/G/A/;
					print CT_OUT $ct_seq, "\n";
					print GA_OUT $ga_seq, "\n";
					$fwd_seq = "NA";
				}
			}
			undef @base;
		}
		close(CH_OUT) if($context);
		if($convert && $pos%150 != 0){
			my ($ct_seq, $ga_seq) = ($fwd_seq, $fwd_seq);
			$ct_seq =~ tr/C/T/;
			$ga_seq =~ tr/G/A/;
			print CT_OUT $ct_seq, "\n";
			print GA_OUT $ga_seq, "\n";
			close(CT_OUT);
			close(GA_OUT);
		}
	}
	close(IN);
}

sub printUsage{
	my $msg = shift;
	print $msg;
	print "Usage 1: genomePrep.pl genome.fa[.gz]\n";
	#print "convert: yes = convert genome, no = don't convert genome\n";
	#print "context: CG = CG only, ALL = all C context\n";
	undef %variant;
	exit 0;
}
main();
