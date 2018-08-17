#!/ifs1/ST_SINGLECELL/USER/jiangrunze/tool/ActivePerl-5.14/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Merge_exp_to_matrix.pl <exp_file.lst> [column_num]
#
#  DESCRIPTION: Merge the gene expression profiles of different samples
#
#       INPUTS: <exp_file.lst> : The list of expression files of all samples
#               [column_num]   : The column index of Reads count, TPM or FPKM
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/04/2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 3;
use File::Basename qw<basename dirname>;

open IN1,shift or die $!;
my $cut = shift;
my $outdir = shift or die $!;
#open OUT1, ">$outdir.ERCC.merge.txt" or die $!;
open OUT2, ">$outdir.gene.merge.txt" or die $!;

$cut -= 1;
my %hash;
my @header;
push @header, "Gene_id";

while (<IN1>){
	chomp;
	open IN2, $_ or die $!;
	my $tmp = <IN2>;
	my $filename = basename $_;
	$filename =~ s/\.genes\.results//;
	push @header, $filename;
	while (<IN2>){
		chomp;
		my @F = split;
		push @{$hash{$F[0]}}, $F[$cut];
	}
}

#print OUT1 join "\t", @header;
print OUT2 join "\t", @header;
#print OUT1 "\n";
print OUT2 "\n";
for (sort keys %hash){
#	if (/^ERCC-/){
#		print OUT1 "$_\t";
#		print OUT1 join "\t", @{$hash{$_}};
#		print OUT1 "\n";
#	}else{
		print OUT2 "$_\t";
		print OUT2 join "\t", @{$hash{$_}};
		print OUT2 "\n";
#	}
}

