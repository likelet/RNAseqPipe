#! /usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Sort_matrix_by_column.pl <Gene_exp.matrix> <Sample_lst>
#
#  DESCRIPTION: Sort matrix by the given sequence of samples
#
#       INPUTS: <Gene_exp.matrix> : Gene expression matrix
#		<Sample_lst>      : Sample annotation file, the first column should be Sample ID
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/07/2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 2;
open IN1,shift or die $!;
open IN2,shift or die $!;
#open OUT,">$ARGV[0]" or die $!;

my %hash;
my %term;;

my $header = <IN1>;
chomp ($header);
my @sample = split /\t/,$header;
for (@sample){
	$term{$_} = ();
}

while (<IN1>){
	chomp;
	my @F = split /\t/;
	for (1..$#F){
		$hash{$F[0]}{$sample[$_]} = $F[$_];
	}
}

my @sort;
while (<IN2>){
	next if /^\t/ or /^Gene/ or /^gene/ or /^Sample_ID/;
	chomp;
	my @F = split /\t/;
	push @sort,$F[0];
}

for(@sort){
	print "\t$_" if exists $term{$_};
}
#print "\t";
#print join "\t", @sort;
print "\n";

for my $gene(sort keys %hash){
	print "$gene";
	for (@sort){
		print "\t$hash{$gene}{$_}" if exists $hash{$gene}{$_};
	}
	print "\n";
}

