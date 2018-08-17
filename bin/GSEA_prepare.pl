#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl GSEA_prepare.pl <Exp.mat> <Sample_info.txt> <Gene_set.gmt> <outdir> [Design.txt]
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
die `pod2text $0` unless @ARGV == 4 or @ARGV == 5;
use strict;
use Math::Combinatorics;

open IN, shift or die $!;
open IN2, shift or die $!;
my $gmt = shift or die $!;
my $outdir = shift or die $!;

my $header = <IN>;
chomp ($header);
my @sample = split /\t/,$header;
my %term;
for (@sample){
        $term{$_} = ();
}

my %hash;
while (<IN>){
        chomp;
        my @F = split /\t/;
        for (1..$#F){
                $hash{$F[0]}{$sample[$_]} = $F[$_];
        }
}

my $header2 = <IN2>;
chomp $header2;
my @G = split /\t/, $header2;
my %index;
for (1..$#G){
	$index{$G[$_]} = $_;
}
my %hash2;
my %hash3;
while(<IN2>){
	chomp;
	my @F = split /\t/;
	push @{$hash2{$F[$index{"Type"}]}}, $F[0];
	$hash3{$F[0]} = $F[$index{"Type"}];
}

my @type = keys %hash2;
my @comb = map{join " ", @$_} combine(2,@type);

sub exp_extact{
	my @tt = @_;
	`mkdir -p $outdir/$tt[0]_vs_$tt[1]`;
	open OUT, "> $outdir/$tt[0]_vs_$tt[1]/exp.mat" or die $!;
        my @sample_list;
	push @sample_list, @{$hash2{$tt[0]}};
	push @sample_list, @{$hash2{$tt[1]}};
        for (@sample_list){
                print OUT "\t$_";
        }
        print OUT "\n";
        for my $gene(sort keys %hash){
	        print OUT "$gene";
	        for (@sample_list){
	                print OUT "\t$hash{$gene}{$_}" if exists $hash{$gene}{$_};
	        }
        	print OUT "\n";
	}
        close OUT;
        open OUT2, "> $outdir/$tt[0]_vs_$tt[1]/anno.cls" or die $!;
        my $count = @sample_list;
        print OUT2 "$count 2 1\n#$tt[0] $tt[1]\n";
        my @line;
        for (@sample_list){
                push @line, $hash3{$_};
        }
        print OUT2 join " ", @line;
        print OUT2 "\n";
        close OUT2;
	print "gseapy gsea -p 20 -d $outdir/$tt[0]_vs_$tt[1]/exp.mat -c $outdir/$tt[0]_vs_$tt[1]/anno.cls -g $gmt -o $outdir/$tt[0]_vs_$tt[1]\n";
#	print "sed \'s/,/\\t/g\' $outdir/$tt[0]_vs_$tt[1]/gseapy.gsea.gene_set.report.csv > $outdir/$tt[0]_vs_$tt[1]/GSEA.gene_set.report.xls\n";
	print "less $outdir/$tt[0]_vs_$tt[1]/gseapy.gsea.gene_set.report.csv | perl -wanle\'\@F=split /,/;my\$a=join \",\", \@F[7..\$#F];print join \"\\t\",\@F[0..6], \$a\' > $outdir/$tt[0]_vs_$tt[1]/GSEA.gene_set.report.xls\n";
}

if ($ARGV[0]){
	open IN, $ARGV[0] or die $!;
	while(<IN>){
		chomp;
		my @rr = split /_vs_/;
		exp_extact(@rr);
	}
}else{
	for my $cc(@comb){
		my @tt = split / /, $cc;
		exp_extact(@tt);
	}
}

