#! /usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Round_off_matrix.pl <Exp_count.mat> > <Exp_count.round_off.mat>
#
#  DESCRIPTION: Round off the non-integral expected_count matrix from RSEM
#  INPUT FILES: <Exp_count.mat> : Merged expression reads count matrix
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/07/2017
#     REVISION: ---
#===============================================================================
=cut
use strict;
open IN, shift or die $!;

my $header = <IN>;
print $header;

while (<IN>){
	chomp;
	my @F = split;
	print $F[0];
	for (1..$#F){
		my $round_off = int($F[$_] + 0.5);
		print "\t$round_off";
	}
	print "\n";
}

