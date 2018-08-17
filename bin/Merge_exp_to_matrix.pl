#!/ifs1/ST_SINGLECELL/USER/jiangrunze/tool/ActivePerl-5.14/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Merge_exp_to_matrix.pl <exp_file.lst> [column_num] > Merged_exp
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
die `pod2text $0` unless @ARGV == 2;
use File::Basename qw<basename dirname>;

open IN1,shift or die $!;
my $cut = shift;
$cut -= 1;
my %hash;
my @header;

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

print "\t";
print join "\t", @header;
print "\n";
for (sort keys %hash){
	print "$_\t";
	print join "\t", @{$hash{$_}};
	print "\n";
}

