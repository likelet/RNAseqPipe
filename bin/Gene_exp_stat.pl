#! /usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Gene_exp_stat.pl <Gene_exp.matrix> > Gene_exp.stat
#
#  DESCRIPTION: Expression statistics of each gene
#
#       INPUTS: <Gene_exp.matrix> : Gene expression matrix
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/07/2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
use strict;
use Statistics::Basic qw(:all);
open IN,shift or die $!;

my $header = <IN>;
chomp ($header);
my @name = split /\t/,$header;
print "$name[0]\tSum\tMean\tMedian\tVariance\tStddev\tCV\t";
print join "\t", @name[1..$#name];
print "\n";
while (<IN>){
	chomp;
	my @F = split;
	my $num = $#F;
	my $gene = shift @F;
	my $mean = mean(@F);
	my $sum = $mean * $num;
	my $median = median(@F);
	my $variance = variance(@F);
	my $stddev = stddev(@F);
	my $CV = "NA";
	$CV = $stddev / $mean unless $mean == 0;
	print "$gene\t$sum\t$mean\t$median\t$variance\t$stddev\t$CV\t";
	print join "\t", @F;
	print "\n";
}

