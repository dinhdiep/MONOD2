#!/usr/bin/perl 
use strict;
use warnings;

my $fqname_file = $ARGV[0];
my $bamfile = $ARGV[1];
my $sample_size = $ARGV[2];

my $max_number = 0;
open(FH, "gunzip -c $fqname_file |") or die("cannot open fqname file\n");
$max_number++ while <FH>;
close(FH);

my %keep;
my %keepfqname;

for(my $i = 0; $i < $sample_size; $i++){
	my $value = int(rand($max_number))+1;
	$keep{$value} = $keep{$value} ? $keep{$value} + 1 : 1;
}

my $increment = 0;
open(FH, "gunzip -c $fqname_file |") or die("cannot open fqname file\n");
while(my $line = <FH>){
	chomp($line);
	$line =~ s/\/1//g;
	$increment++;
	if($keep{$increment}){
		$keepfqname{$line}=$keep{$increment};
		#print $line, "\n";
	}
}
close(FH);

open(FH, "samtools view $bamfile |") or die("cannot convert bamfile\n");
while(my $line = <FH>){
	chomp($line);
	my @tmp = split "\t", $line;
	next if(!$keepfqname{$tmp[0]});
	for(my $i = 0; $i < $keepfqname{$tmp[0]}; $i++){
		print $line, "\n";
	}
}
close(FH);

