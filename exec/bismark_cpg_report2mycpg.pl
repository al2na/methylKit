#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $inputfile;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
use File::Basename;
my $tag = 'myCpG';
my $add_chr;
my $context = 'CpG';

GetOptions(

           'i|input|inputfile:s'  => \$inputfile,
           'tag:s'                => \$tag,
           'context:s'            => \$context,
           'add_chr'              => \$add_chr,
           'debug'                => \$debug,
           'verbose'              => \$verbose,
           'simulate'             => \$simulate,

          );

my @suffixlist = '.CpG_report.txt.gz';
open IN, "gunzip -c $inputfile |" or die $!;
my ($name,$path,$suffix) = fileparse($inputfile,@suffixlist);
my $outfile = $path . $name . ".$tag.txt";
open OUT, ">$outfile" or die $!;
# header
print OUT "chrBase" . "\t"
            . "chr" . "\t"
            . "base" . "\t"
            . "strand" . "\t"
            . "coverage" . "\t"
            . "freqC" . "\t"
            . "freqT" . "\n";

while (<IN>) {
  chomp $_;
  my ($chr,$base,$sign_strand,$met_count,$unm_count,$c_context,$trinuc_context) = split("\t",$_);
  next if (0 == $met_count && 0 == $unm_count);
  if ($add_chr) { my $tmp = $chr; $chr = 'chr' . $tmp; }
  my $chrBase = "$chr.$base";
  my $strand = 'F'; $strand = 'R' if ($sign_strand eq '-');
  my $coverage = $met_count + $unm_count;
  my $freqC = "0.0"; $freqC = sprintf("%.02f", (100*$met_count)/$coverage) if ($met_count > 0);
  my $freqT = "0.0"; $freqT = sprintf("%.02f", (100*$unm_count)/$coverage) if ($unm_count > 0);
  print OUT "$chrBase"  . "\t"
    . "$chr" . "\t"
    . "$base" . "\t"
    . "$strand" . "\t"
    . "$coverage" . "\t"
    . "$freqC" . "\t"
    . "$freqT" . "\n";
}
close IN;
close OUT;

print "outfile:\n";
print "$outfile\n";

1;

$DB::single=1;1;
$DB::single=1;1;

# bismark_cpg_report2mycpg.pl
#
# Cared for by Albert Vilella <avilella@gmail.com>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

bismark_cpg_report2mycpg.pl - DESCRIPTION 

=head1 SYNOPSIS

Transforms from .CpG_report.txt.gz to .myCpG.txt format as below:

The genome-wide cytosine report(optional) is tab-delimited in the following format:

(1-basedcoords):<chromosome>  <position>  <strand>  <count methylated>  <count unmethylated>  <C-context> <trinucleotide context>

The .myCpG.txt format is as follows:

# chrBase	chr	base	strand	coverage	freqC	freqT
# chr21.9826907	chr21	9826907	F	96	18.75	81.25
# chr21.9853326	chr21	9853326	F	16	87.50	12.50
# chr21.9853296	chr21	9853296	F	18	88.89	11.11
# chr21.9860126	chr21	9860126	F	83	100.00	0.00
# chr21.9906663	chr21	9906663	R	14	92.86	7.14
# chr21.9906677	chr21	9906677	R	14	78.57	21.43
# chr21.9906700	chr21	9906700	F	30	26.67	73.33
# chr21.9906954	chr21	9906954	F	30	86.67	13.33
# chr21.9906644	chr21	9906644	F	23	95.65	4.35
# chr21.9906873	chr21	9906873	F	15	26.67	73.33
# chr21.9906634	chr21	9906634	F	23	69.57	30.43
# chr21.9906704	chr21	9906704	F	18	94.44	5.56
# chr21.9906655	chr21	9906655	R	14	92.86	7.14
# chr21.9906645	chr21	9906645	R	13	100.00	0.00

=head1 DESCRIPTION

GetOptions(

           'i|input|inputfile:s'  => \$inputfile,
           'tag:s'                => \$tag,
           'debug'                => \$debug,
           'verbose'              => \$verbose,
           'simulate'             => \$simulate,

=head1 AUTHOR - Albert Vilella

Email avilella@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

