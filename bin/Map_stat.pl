#!/usr/bin/perl -w
die "perl $0 <Base_dir>" unless @ARGV == 1;
use File::Basename; 

my $outdir = shift or die $!;
my @file = `ls $outdir/2.ReadMap/*/Log.final.out`;
print "Sample_ID\tTotal_reads_count\tUniq_map_reads\tUniq_map_rate\tMuti_map_reads\tMuti_map_count\n";
for my $f(@file){
	$f =~ /ReadMap\/(\S+)\//;
	my $name = $1;
	open IN, $f or die $!;
	print $name;
	while (<IN>){
		chomp;
		my @F = split /\t/;
		print "\t$F[1]" if /Number of input reads/ or /Uniquely mapped reads number/ or /Uniquely mapped reads %/ or /Number of reads mapped to multiple loci/ or /\% of reads mapped to multiple loci/;
	}
	print "\n";
}

