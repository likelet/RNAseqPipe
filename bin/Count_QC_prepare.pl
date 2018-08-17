#!/ifs1/ST_SINGLECELL/USER/jiangrunze/tool/ActivePerl-5.14/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Count_QC_prepare.pl <Out_dir>
#
#  DESCRIPTION: Preparation of qualimap-CountQC
#
#       INPUTS: <Out_dir> : The root analysis derectory of this project
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/09/2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
use File::Basename qw<basename dirname>;

my $outdir = shift or die $!;
my @file = `ls $outdir/3.Expression`;
my %hash;
my @header;

for my $f(@file){
	chomp $f;
	open IN2, "$outdir/3.Expression/$f/$f.genes.results" or die $!;
	my $tmp = <IN2>;
	$f =~ s/\.genes\.results//;
	push @header, $f;
	while (<IN2>){
		chomp;
		my @F = split;
		push @{$hash{$F[0]}}, $F[4];
	}
}

open OUT1, ">$outdir/tmp/List_for_count_QC.txt" or die $!;
open OUT2, ">$outdir/tmp/Count_matrix_for_count_QC.txt" or die $!;
my $n = 2;
for (@header){
	chomp;
	print OUT1 "$_\t1\t$outdir/tmp/Count_matrix_for_count_QC.txt\t$n\n";
	$n ++;
}

for (sort keys %hash){
	print OUT2 "$_\t";
	print OUT2 join "\t", @{$hash{$_}};
	print OUT2 "\n";
}

