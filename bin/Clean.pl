#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Clean.pl <project_path> 
#
#  DESCRIPTION: 
#
#  INPUT FILES: 
#
# REQUIREMENTS:
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: //2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
use strict;

my $outdir = shift or die $!;
# 2.ReadMap
`rm -rf $outdir/2.ReadMap/*/Log.* $outdir/2.ReadMap/*/_* $outdir/2.ReadMap/*/SJ.out.tab`;
# 3.Expression
`rm -rf $outdir/3.Expression/*/*.stat $outdir/3.Expression/*/*.time`;
# 4.MapQC

# 5.Fusion
`rm -rf $outdir/5.Fusion/*/_* $outdir/5.Fusion/*/star-fusion.preliminary`;
# 6.SNV
`rm -rf $outdir/6.SNV/*/Variant.filt.avinput $outdir/6.SNV/*/Variant.filtered.avinput $outdir/6.SNV/*/Variant.filt.vcf.idx $outdir/6.SNV/*/Variant.vcf.idx`;
# 7.Further_analysis
`rm -rf $outdir/7.Further_analysis/4.GSEA/*/anno.cls $outdir/7.Further_analysis/4.GSEA/*/exp.mat $outdir/7.Further_analysis/4.GSEA/*/gseapy.gsea.gene_set.report.csv $outdir/7.Further_analysis/4.GSEA/*/gseapy.gsea.go_kegg_reactome.v6.1.symbols.log`;


