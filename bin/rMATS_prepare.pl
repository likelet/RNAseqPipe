#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl rMATS_prepare.pl <Sample_info.txt> <outdir> [reverse]
#
#  DESCRIPTION: GSEA preparision for paired-with sample types
#
#  INPUT FILES: 
#
# REQUIREMENTS:
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 4/2/2018
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 2 or @ARGV == 3;
use strict;
use Math::Combinatorics;

open IN, shift or die $!;
my $outdir = shift or die $!;;

my $header = <IN>;
my %hash;
my %hash2;
while(<IN>){
	chomp;
	my @F = split /\t/;
	push @{$hash{$F[1]}}, $F[0];
	$hash2{$F[0]} = $F[1];
}

my @type = keys %hash;
my %bam_path;
for my $t(@type){
	open OUT, ">$outdir/tmp/$t\_bam.txt" or die $!;
	my @bam;
	for my $s(@{$hash{$t}}){
		my $bam_file = "$outdir/2.ReadMap/$s/Aligned.sortedByCoord.out.bam";
		push @bam, $bam_file;
	}
	print OUT join ",", @bam;
	$bam_path{$t} = join ",", @bam;
	print OUT "\n";
	close OUT;
}

my @comb = map{join " ", @$_} combine(2,@type);
open OUT2, ">$outdir/Shell/Conbine_analysis/5.Alt_splicing.sh" or die $!;
for my $cc(@comb){
	my @tt = split / /, $cc;
	`mkdir -p $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]`;
	print OUT2 "# Differential alt-splicing detection\n";
	print OUT2 "python /disk/soft/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 $outdir/tmp/$tt[1]_bam.txt --b2 $outdir/tmp/$tt[0]_bam.txt --gtf /disk/database/human/hg38/Gencode/GRCh38_gencode_v24_CTAT_lib_Mar292017/GRCh38.gtf --od $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0] -t paired --readLength 151 --cstat 0.0001 --nthread 10\n";
	print OUT2 "# Alt-splicing filtering\n";
	print OUT2 "awk \'\$20 < 0.05\' $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/A3SS.MATS.JCEC.txt > $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/A3SS.MATS.JCEC.filt.txt\n";
	print OUT2 "awk \'\$20 < 0.05\' $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/RI.MATS.JCEC.txt > $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/RI.MATS.JCEC.filt.txt\n";
	print OUT2 "awk \'\$20 < 0.05\' $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/A5SS.MATS.JCEC.txt > $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/A5SS.MATS.JCEC.filt.txt\n";
	print OUT2 "awk \'\$20 < 0.05\' $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/MXE.MATS.JCEC.txt > $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/MXE.MATS.JCEC.filt.txt\n";
	print OUT2 "awk \'\$20 < 0.05\' $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/SE.MATS.JCEC.txt > $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/SE.MATS.JCEC.filt.txt\n";
	print OUT2 "# Ploting alt-splicing\n";
	print OUT2 "rmats2sashimiplot --b1 $bam_path{$tt[1]} --b2 $bam_path{$tt[0]} -e $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/A3SS.MATS.JCEC.filt.txt -t A3SS --l1 $tt[1] --l2 $tt[0] -o $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/Figure_A3SS &\n";
	print OUT2 "rmats2sashimiplot --b1 $bam_path{$tt[1]} --b2 $bam_path{$tt[0]} -e $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/RI.MATS.JCEC.filt.txt -t RI --l1 $tt[1] --l2 $tt[0] -o $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/Figure_RI &\n";
	print OUT2 "rmats2sashimiplot --b1 $bam_path{$tt[1]} --b2 $bam_path{$tt[0]} -e $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/A5SS.MATS.JCEC.filt.txt -t A5SS --l1 $tt[1] --l2 $tt[0] -o $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/Figure_A5SS &\n";
	print OUT2 "rmats2sashimiplot --b1 $bam_path{$tt[1]} --b2 $bam_path{$tt[0]} -e $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/MXE.MATS.JCEC.filt.txt -t MXE --l1 $tt[1] --l2 $tt[0] -o $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/Figure_MXE &\n";
	print OUT2 "rmats2sashimiplot --b1 $bam_path{$tt[1]} --b2 $bam_path{$tt[0]} -e $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/SE.MATS.JCEC.filt.txt -t MATS --l1 $tt[1] --l2 $tt[0] -o $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]/Figure_MATS &\n";
	if ($ARGV[0]){
#		`mkdir -p $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0]`;
#		print OUT2 "python /disk/soft/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 $outdir/tmp/$tt[1]_bam.txt --b2 $outdir/tmp/$tt[0]_bam.txt --gtf /disk/database/human/hg38/Gencode/GRCh38_gencode_v24_CTAT_lib_Mar292017/GRCh38.gtf --od $outdir/7.Further_analysis/5.Alt_splicing/$tt[1]_vs_$tt[0] -t paired --readLength 151 --cstat 0.0001 --nthread 10\n";
	}
}	
