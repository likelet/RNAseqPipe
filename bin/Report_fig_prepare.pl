#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Report_fig_prepare.pl <outdir>
#
#  DESCRIPTION: Generate png figures for report.
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

`mkdir -p $outdir/tmp/Images`;

# Replot fusion circos plot
`grep Rscript $outdir/Shell/*/5.fusion.sh | perl -wanle '\@F=split /:/;print \$F[1]' | xargs -iCMD -P10 bash -c CMD`;

# Mapping stat
if (-e "$outdir/7.Further_analysis/1.QC/MapQC/Map_stat.pdf"){
	`convert $outdir/7.Further_analysis/1.QC/MapQC/Map_stat.pdf $outdir/tmp/Images/Map_stat.png`;
}else{
	print "$outdir/7.Further_analysis/1.QC/MapQC/Map_stat.pdf is not exist\n";
}
# Expression
if (-e "$outdir/7.Further_analysis/2.Exp_Merge/Exp_TPM.all_0.TPM_dis.pdf"){
	`convert $outdir/7.Further_analysis/2.Exp_Merge/Exp_TPM.all_0.TPM_dis.pdf $outdir/tmp/Images/Exp_TPM.all_0.TPM_dis.png`;
}else{
	print "$outdir/7.Further_analysis/2.Exp_Merge/Exp_TPM.all_0.TPM_dis.pdf is not exist\n";
}
# DEG
open GP, "$outdir/tmp/Design.txt" or die $!;
while(<GP>){
	chomp;
	my $group = $_;
	if (-e "$outdir/7.Further_analysis/3.DEG/$group.deseq.Plot.pdf"){
		`convert $outdir/7.Further_analysis/3.DEG/$group.deseq.Plot.pdf $outdir/tmp/Images/$group.deseq.Plot.png`;
	}else{
		print "$outdir/7.Further_analysis/3.DEG/$group.deseq.Plot.pdf is not exist\n";
	}
	if(-e "$outdir/tmp/image/$group.deseq.Plot.png"){
		`mv $outdir/tmp/image/$group.deseq.Plot.png $outdir/tmp/Images/$group.deseq.Plot-0.png`;
	}
#Function enrich
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_BP_enrich_plot1.pdf"){
		`convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_BP_enrich_plot1.pdf $outdir/tmp/Images/$group.up.GO_BP.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_BP_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_BP_enrich_plot1.pdf"){
		`convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_BP_enrich_plot1.pdf $outdir/tmp/Images/$group.down.GO_BP.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_BP_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_MF_enrich_plot1.pdf"){
                `convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_MF_enrich_plot1.pdf $outdir/tmp/Images/$group.up.GO_MF.png`;
        }else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_MF_enrich_plot1.pdf\n";
	}       
        if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_MF_enrich_plot1.pdf"){
                `convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_MF_enrich_plot1.pdf $outdir/tmp/Images/$group.down.GO_MF.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_MF_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_CC_enrich_plot1.pdf"){
                `convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_CC_enrich_plot1.pdf $outdir/tmp/Images/$group.up.GO_CC.png`;
        }else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.GO_CC_enrich_plot1.pdf\n";
	}       
        if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_CC_enrich_plot1.pdf"){
                `convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_CC_enrich_plot1.pdf $outdir/tmp/Images/$group.down.GO_CC.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.GO_CC_enrich_plot1.pdf\n";
        }
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.KEGG_enrich_plot1.pdf"){
		`convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.KEGG_enrich_plot1.pdf $outdir/tmp/Images/$group.up.KEGG.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.KEGG_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.KEGG_enrich_plot1.pdf"){
                `convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.KEGG_enrich_plot1.pdf $outdir/tmp/Images/$group.down.KEGG.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.KEGG_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.Reactome_enrich_plot1.pdf"){
		`convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.Reactome_enrich_plot1.pdf $outdir/tmp/Images/$group.up.Reactome.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.Reactome_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.Reactome_enrich_plot1.pdf"){
                `convert $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.Reactome_enrich_plot1.pdf $outdir/tmp/Images/$group.down.Reactome.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.Reactome_enrich_plot1.pdf\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.KEGG_enrich.xls"){
		open IN, "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.KEGG_enrich.xls";
		my $header_up = <IN>;
		my $KEGG_up = <IN>;
		my $KEGG_up_1 = (split /\t/, $KEGG_up)[0];
		`cp $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$KEGG_up_1.$group.deseq.up_regulate_enrich_all.png $outdir/tmp/Images/$group.up.KEGG.path.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.up_regulate.KEGG_enrich.xls\n";
	}
	if (-e "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.KEGG_enrich.xls"){
                open IN, "$outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.KEGG_enrich.xls";
                my $header_down = <IN>;
                my $KEGG_down = <IN>;
                my $KEGG_down_1 = (split /\t/, $KEGG_down)[0];
                `cp $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$KEGG_down_1.$group.deseq.down_regulate_enrich_all.png $outdir/tmp/Images/$group.down.KEGG.path.png`;
	}else{
		print "No $outdir/7.Further_analysis/3.DEG/Function_enrich/$group/$group.deseq.down_regulate.KEGG_enrich.xls\n";
	}
# GSEA
	if (-e "$outdir/7.Further_analysis/4.GSEA/$group/GSEA.gene_set.report.xls"){
		open IN, "$outdir/7.Further_analysis/4.GSEA/$group/GSEA.gene_set.report.xls";
		my $header_gsea = <IN>;
		my $GSEA = <IN>;
		my @gsea = split /\t/, $GSEA;
		if ($gsea[4]<0.05){
			`convert $outdir/7.Further_analysis/4.GSEA/$group/$gsea[0].gsea.pdf $outdir/tmp/Images/$group.GSEA.png`;
		}
	}else{
		print "No $outdir/7.Further_analysis/4.GSEA/$group/GSEA.gene_set.report.xls\n";
	}
}
	
my @sample = `ls $outdir/5.Fusion`;
for my $ss(@sample){
	chomp $ss;
	`convert $outdir/5.Fusion/$ss/Fusion_CircosPlot.pdf $outdir/tmp/Images/$ss.fusion_circos.png`;
=t
	if (-e "$outdir/5.Fusion/$ss/star-fusion.fusion_predictions.abridged.tsv"){
		my $nrow = `wc -l $outdir/5.Fusion/$ss/star-fusion.fusion_predictions.abridged.tsv`;
		print "No results in $outdir/5.Fusion/$ss/star-fusion.fusion_predictions.abridged.tsv\n" and next if $nrow =~ /^1 /;
		`mv $outdir/5.Fusion/$ss/star-fusion.fusion_predictions.abridged.tsv $outdir/tmp/Images/$ss.fusion.txt`;
	}else{
		print "No $outdir/5.Fusion/$ss/star-fusion.fusion_predictions.abridged.tsv\n";
	}
=cut
}	
