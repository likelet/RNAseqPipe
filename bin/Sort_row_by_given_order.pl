#! /usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Specific_gene_extract.pl <Gene_exp.matrix> <Gene_lst>
#
#  DESCRIPTION: Expression statistics of each gene
#
#       INPUTS: <Gene_exp.matrix> : Gene expression matrix
		<Gene_lst>        : Gene list you want to extract
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/07/2017
#     REVISION: ---
#===============================================================================
=cut
#die `pod2text $0` unless(@ARGV == 4 or @ARGV == 5);
open IN1,shift or die $!;
open IN2,shift or die $!;
#open OUT,">$ARGV[0]" or die $!;
my %genes;
my %hash;

my $head=<IN1>;
print $head;
while(<IN1>){
        my @F=split /\t/;
	$genes{$F[0]} = $_;
}

while(<IN2>){
	chomp;
	my @F = split;
	print $genes{$F[0]} if exists $genes{$F[0]};
}

