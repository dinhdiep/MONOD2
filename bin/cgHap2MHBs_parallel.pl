#!/usr/bin/perl -w
# cgHap2MHBs_parallel.pl

use Parallel::Simple qw(prun);

if(!$ARGV[3]){
  print "Usage:\n";
  print " perl $0 [hapInfo file] [min R2] [snpBed file] [outPrefix]\n\n";
  exit(1);
}

my $hapInfo=$ARGV[0];
my $minR2=$ARGV[1];
my $snpBedFile=$ARGV[2];
my $outNamePrefix=$ARGV[3];

my @args = ("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu");

my $num_processes = scalar(@args);

sub run_process{
#  my $hap_info_file = 
#  my $bin_file = 
#  my $min_r2 = 
#  my $cg_snp_bed_file =
  my $file = shift;
  my $cmd = "bin/hapInfo_maskSNPs2mld_block.pl hapinfo_$file.tmp $minR2 $snpBedFile > $outNamePrefix.$file";
  system($cmd);
  unlink("hapinfo_$file.tmp");
}

my %fileHandles;
# split haps into 12 files
for(my $i = 0; $i < scalar(@args); $i++){
  my $fname = "hapinfo_" . $args[$i] . ".tmp";
  open( $fileHandles{$args[$i]} , ">$fname");
}

my ($cur_index, $prev_index) = ("NA", "NA");
my $cur_i_index = 0;
open(INFILE, "$hapInfo") || die("Error reading $hapInfo\n");
while(my $line = <INFILE>){
  my @temp = split "\t", $line;
  my $cur_index = $temp[0];
  if($prev_index ne "NA" and $cur_index eq $prev_index){
     print { $fileHandles{ $args[$cur_i_index % $num_processes] } } $line if(length($temp[1]) >= 3);
  }else{
     $cur_i_index++;
     print { $fileHandles{ $args[$cur_i_index % $num_processes] } } $line if(length($temp[1]) >= 3);
  }
  $prev_index = $cur_index;
}
close(INFILE);

for(my $i = 0; $i < scalar(@args); $i++){
  close( $fileHandles{$args[$i]} );
}

prun( 
   [ \&run_process, $args[0] ],
   [ \&run_process, $args[1] ],
   [ \&run_process, $args[2] ],
   [ \&run_process, $args[3] ],
   [ \&run_process, $args[4] ],
   [ \&run_process, $args[5] ],
   [ \&run_process, $args[6] ],
   [ \&run_process, $args[7] ], 
   [ \&run_process, $args[8] ], 
   [ \&run_process, $args[9] ], 
   [ \&run_process, $args[10] ], 
   [ \&run_process, $args[11] ], 
) or die( Parallel::Simple::errplus() );

print "all done!\n";

