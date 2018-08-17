#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Z_score.pl <Gene_exp.matrix> > Gene_exp.zscore.matrix
#
#  DESCRIPTION: Convert expression matrix into zscore matrix
#
#  INPUT FILES: <Gene_exp.matrix> : Gene expression matrix
#
# REQUIREMENTS: Perl mudule Statistics::Basic
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: //2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
use strict;
use Statistics::Basic qw(:all);
open IN,shift or die $!;

my $header = <IN>;
print $header;

while(<IN>){
	chomp;
	my @F = split;
	my $mean = mean(@F[1..$#F]);
	my $stddev = stddev(@F[1..$#F]);
	for my $i(1..$#F){
		$F[$i] = ($F[$i] - $mean) / $stddev;
	}
	print join "\t", @F, "\n";
}
