#! /usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Ensem2Symbol2.pl  <exp_matrix_from_DAtools>
#
#  DESCRIPTION: Convert Ensembel_ID into Gene_symbol, collapse count by mean
#
#       INPUTS: <anno_file>  : Ensembol_id to GeneSymbol_id file
#               <exp_matrix> : Expression matrix file
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/04/2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;

my $input = shift or die $!;

open IN1, $input or die $!;


my %hash;
my %mark;

while (<IN1>){
	chomp;
	my @F = split /\t/;
	$hash{$F[0]} = $F[1];
}

open IN2, $input or die $!;

my $header = <IN2>;
my @G = split /\t/, $header;
print join "\t", @G[1,3..$#G];
while (<IN2>){
	my @F = split;
	next unless exists $hash{$F[0]};
	if (exists $mark{$hash{$F[0]}}){
		for (3..$#F){
			$hash{$F[0]}[$_] += $F[$_];
			$hash{$F[0]}[$_] = "0.00" if $hash{$F[0]}[$_] == 0;
		}
	}else{
		$mark{$hash{$F[0]}} = ();
		@{$hash{$F[0]}} = @F;
		${$hash{$F[0]}}[0] = $hash{$F[0]};
	}
}
for (sort keys %mark){
	print join "\t", @{$_}[1,3..$#{$_}];
	print "\n";
}
