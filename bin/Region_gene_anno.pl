#!/usr/bin/perl -w
die "perl $0 <CNV.bed> <Gene.bed>" unless @ARGV == 2;
open CNV, shift or die $!;
open ANNO,shift or die $!;
#open OUT,">$ARGV[0]" or die $!;

my %hash;

my $tmp_chr = 1;
my $tmp_start = 1;
my $tmp_end = 1;
my $tmp_gene = "";

my $header = <CNV>;
chomp $header;
print "$header\tGenes\tGene_num\n";

while (<CNV>){
	chomp;
	next if /^chrX/ or /^chrY/ or /^chrM/ or /^GL/;
	my @F = split;
	$F[0] =~ /chr(\d+)/;
	my $chr_cnv = $1;
	print $_."\t";
	my $count = 0;
	if (($chr_cnv == $tmp_chr and $F[1] > $tmp_start) or $chr_cnv > $tmp_chr){
	}elsif(($chr_cnv == $tmp_chr and $F[2] < $F[1]) or $chr_cnv < $tmp_chr){
		print "\t0\n";
		next;
	}else{
		print "$tmp_gene;";
		$count ++;
	}
	while (<ANNO>){
		chomp;
		next if /^chrX/ or /^chrY/ or /^chrM/ or /^GL/;
		my @G = split;
		$G[0] =~ /chr(\d+)/;
		my $chr_gene = $1;
		if (($chr_cnv == $chr_gene and $F[1] > $G[2]) or $chr_cnv > $chr_gene){
			next;
		}elsif(($chr_cnv == $chr_gene and $F[2] < $G[1]) or $chr_cnv < $chr_gene){
			$tmp_chr = $chr_gene;
			$tmp_start = $G[1];
			$tmp_end = $G[2];
			$tmp_gene = $G[3];
			last;
		}else{
			print "$G[3];";
			$count ++;
		}
	}
	print "\t$count\n";
}

