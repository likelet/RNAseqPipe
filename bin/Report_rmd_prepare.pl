#!perl -w
=head1 #===============================================================================
#        USAGE: perl Report_rmd_prepare.pl <outdir> > report.rmd
#
#  DESCRIPTION: Gernerate report.rmd for whole project.
#
#       INPUTS: 
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 05/02/2018
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;

my $outdir = shift;
my $pipedir = "/disk/zxl/tools/RNA-seq-Pipe";

print "---
title: \"转录组生物信息分析报告\"
author: \"张小龙\"
date: \"`r format(Sys.time(), '%Y年%m月%d日')`\"
output:
  html_document: 
    toc: yes
    toc_depth: 2
    toc_float:
      collapsed: false
      smooth_scroll: false
fontsize: 20pt
params:
  outdir: \".\"
---

<style type=\"text/css\">

h1 { /* Header 1 */
 font-size: 28px;
 color: DarkBlue;
}
h2 { /* Header 2 */
 font-size: 22px;
 color: DarkBlue;
}
h3 { /* Header 3 */
 font-size: 20px;
 color: DarkBlue;
}
h4 { /* Header 4 */
 font-size: 20px;
 color: DarkBlue;
}
h5 { /* Header 5 */
 font-size: 20px;
 color: DarkBlue;
}
h6 { /* Header 6 */
 font-size: 20px;
 color: DarkBlue;
}
h7 { /* Header 7 */
 font-size: 20px;
 color: DarkBlue;
}
</style>

```{r setup, include=FALSE}
knitr::opts_chunk\$set(echo = TRUE, fig.align=\"center\")
```

<img src=\"/disk/zxl/tools/RNA-seq-Pipe/Images/SYSUCC_logo.png\" 
style=\"position:absolute;top:10px;right:30px;width:3.5cm\" />

# 一、建库测序流程

从RNA样品提取到最终数据获得，样品检测、建库、测序等每一环节都会直接影响数据的数量和质量，从而影响后续信息分析的结果。建库测序的流程图如下：

```{r out.width = \"25%\", echo=FALSE}
img1_path = \"$pipedir/Images/sequencing_flow.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图1.1 建库测序流程</center>


## 1.	RNA提取与检测

采用标准提取方法从组织或细胞中提取RNA，随后对RNA样品进行严格质控，质控标准主要包括以下4个方面：

- 琼脂糖凝胶电泳：分析样品RNA完整性及是否存在DNA污染；
- NanoPhotometer spectrophotometer：检测RNA纯度(OD260/280及OD260/230比值)；
- Qubit2.0 Fluorometer：RNA浓度精确定量；
- Agilent 2100 bioanalyzer：精确检测RNA完整性。

## 2.	文库构建与质检

### 2.1文库构建

mRNA的获取主要有两种方式：一是利用真核生物大部分mRNA都带有polyA尾的结构特征，通过Oligo(dT)磁珠富集带有polyA尾的mRNA。二是从总RNA中去除核糖体RNA，从而得到mRNA。随后在NEB Fragmentation Buffer中用二价阳离子将得到的mRNA随机打断，按照NEB普通建库方式或链特异性建库方式进行建库。

**NEB普通建库：**以片段化的mRNA为模版，随机寡核苷酸为引物，在M-MuLV逆转录酶体系中合成cDNA第一条链，随后用RNaseH降解RNA链，并在DNA polymerase I 体系下,以dNTPs为原料合成cDNA第二条链。纯化后的双链cDNA经过末端修复、加A尾并连接测序接头(1)，用AMPure XP beads筛选200bp左右的cDNA， 进行PCR扩增并再次使用AMPure XP beads纯化PCR产物，最终获得文库。建库原理如下图左所示。

**链特异性建库：**逆转录合成cDNA第一条链方法与NEB普通建库方法相同，不同之处在于合成第二条链时，dNTPs中的dTTP由dUTP取代，之后同样进行cDNA末端修复、加A尾、连接测序接头和长度筛选，然后先使用USER酶降解含U的cDNA第二链再进行PCR扩增并获得文库。链特异性文库具有诸多优势，如相同数据量下可获取更多有效信息；能获得更精准的基因定量、定位与注释信息；能提供反义转录本及每一isoform中单一exon的表达水平。建库原理如下图右所示。

```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$pipedir/Images/library_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图1.2 建库原理</center>


**注：**测序接头：包括P5/P7，index和Rd1/Rd2 SP三个部分（如下图所示）。其中P5/P7是PCR扩增引物及flow cell上引物结合的部分，index提供区分不同文库信 息的Rd1/Rd2，SP即read1/read2 sequence primer，是测序引物结合区域，测序过程理论上由Rd1/Rd2 SP向后开始进行。


### 2.2文库质检

文库构建完成后，先使用Qubit2.0 Fluorometer进行初步定量，稀释文库至1.5ng/ul，随后使用Agilent 2100 bioanalyzer对文库的insert size进行检测，insert size符合预期后，qRT-PCR对文库有效浓度进行准确定量(文库有效浓度高于2nM)，以保证文库质量。

## 3.	上机测序

库检合格后，把不同文库按照有效浓度及目标下机数据量的需求pooling后进行Illumina HiSeq测序。测序的基本原理是边合成边测序（Sequencing by Synthesis)。在 测序的flow cell中加入四种荧光标记的dNTP、DNA聚合酶以及接头引物进行扩增，在每一个测序簇延伸互补链时，每加入一个被荧光标记的dNTP就能释放出相对应的荧光，测序仪通过捕获荧光信号，并通过计算机软件将光信号转化为测序峰，从而获得待测片段的序列信息。测序过程如下图所示。


```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$pipedir/Images/sequencing_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图1.3 Illumina测序过程</center>


# 二、信息分析流程

## 1. 标准分析流程

RNA-seq的核心是基因表达差异的显著性分析，使用统计学方法，比较两个条件或多个条件下的基因表达差异，从中找出与条件相关的特异性基因，然后进一步分析这些特异性基因的生物学意义，分析过程包括质控，比对，定量，差异显著性分析，功能富集六个环节，如下图所示。另外可变剪接，变异位点，融合基因也是RNA-seq的重要分析内容。

```{r out.width = \"50%\", echo=FALSE}
img1_path = \"$pipedir/Images/basic_analysis_pipe.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图2.1 基本流程分析内容</center>


## 2. 高级定制化分析

在基本流程分析的基础上，通过样本表达谱对样本进行聚类与分型，同时整合TCGA数据，扩充样本量，使得结果更加可靠。通过临床关联及生存分析，寻找与临床表型潜在关联的基因标记，可用于进一步的功能实验。通过基因共表达分析，对基因进行模块化区分，并通过与临床信息的关联寻找具有生物学意义的共表达模块，并通过网络的构建，寻找模块中核心基因（Hub-gene）。

```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$pipedir/Images/personalized_analysis_pipe.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图2.2 高级定制化分析内容</center>


# 三、结果展示

## 1.原始数据质控分析

### 1.1 原始数据说明

高通量测序（如Illumina HiSeq PE125/PE150）下机得到的原始图像文件经CASAVA碱基识别转化为测序读段（Sequenced Reads），以FASTQ格式存储。FASTQ是一种存储生物序列及相应质量值的常用文本格式，每条测序读段由四行组成，如下图所示。

```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$pipedir/Images/fastq_file.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```
<center>图3.1.1 fastq格式文件内容</center>


第一行为测序读段的标识信息，通常以“@”开头，随后为Illumina测序标识别符和描述信息。第二行为测序读段的碱基序列。第三行通常以“+”开头，存储与第一行相同的信息或为缺省值。第四行为对应碱基的测序质量值，该行字符为第二行对应碱基的质量值加上33后转换为ASCII码，逆向转化即可直观得到每个碱基的质量信息。

### 1.2 测序质量分布

RNA-seq技术测序错误率分布存在如下特点：1.测序读段前几个碱基的错误率较高，这一般是由于测序仪在测序前期稳定性较差，图像识别质量低，另外接头空载，也会引起错误率升高。2.随着测序读段长度的延伸，测序错误率呈上升趋势，这可能是由于测序过程中，每个cycle在荧光基团淬灭，去3'端保护基团时，没能完全去除，导致在延伸时滞留，或者是加入了无3'端保护的碱基，导致延伸超前，滞留或超前引起延伸步调不一致。这是一个累积的过程，越到后面，滞留或超前的累积越多，测序错误率也就越高。各样本测序质量值分布如图1.1所示。绿色区域为理想质量分布区域。

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_fastqc_per_base_sequence_quality_plot_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.1.2 样本测序质量值分布</center>


### 1.3 GC含量分布

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_fastqc_per_sequence_gc_content_plot_Percentages.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.1.3 样本测序GC含量分布</center>


### 1.4 重复序列水平检测

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_fastqc_sequence_duplication_levels_plot_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.1.4 测序重复水平分布</center>


### 1.5 接头污染水平检测

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_fastqc_adapter_content_plot_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.1.5 测序接头污染情况分布</center>


## 2.比对分析
将clean reads比对到基因组或转录组是后续分析的基础。对于基因注释信息比较完善的物种如人和小鼠），如果只研究表达差异显著性基因，可直接将reads比对到转录组。如果研究可变剪接，变异，融合基因等内容，需在基因组水平进行比对。若选择比对到基因组，所选用的比对软件要求能够进行splice，否则会丢失很多有效的junction reads，这是因为可变剪接和polyA尾在真核生物中普遍存在，跨越可变剪接位点或ployA尾边缘的reads无法完全比对到基因组上。本分析流程采用STAR软件对RNA-seq测序数据进行比对分析，STAR采用Maximal Mappable Prefix（MMP）搜索方法，可以对junction reads进行精确定位，如下图所示，其综合性能在同类比对软件中表现较为突出。

```{r out.width = \"50%\", echo=FALSE}
img1_path = \"$pipedir/Images/star_align_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.2.1 STAR对junction reads进行精确定位</center>


### 2.1 比对结果说明

比对输出结果一般为标准的SAM文件，由sanger制定，以TAB为分隔符的一种文本格式，主要用于测序序列比对到基因组上的结果展示。BAM文件是SAM文件的二进制格式，因其占据存储资源小，运算速度快，很多比对软件以BAM格式作为标准输出，但想要查看其文件内容，需要借助samtools工具。SAM文件一般包括两个部分，注释信息（header section）和比对结果（alignment section）。注释信息可有可无，都是以“\@”开头，用不同的tag表示不同的信息，主要有“\@HD”（说明符合标准的版本、对比序列的排列顺序）；“\@SQ”（参考序列说明）；“\@RG”（比对上的reads说明）；\@PG(使用的程序说明)；\@CO(其它任意的说明信息)。比对结果部分，每一行表示一个片段（segment）的比对信息，包括11个必须的字段（mandatory fields）和部分可选的字段，字段之间用TAB键分隔。字段的顺序是固定的，不可用时，根据字段定义，可以用“0”或者“*”代替。STAR软件的输出结果如下表所示：

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$pipedir/Images/Bam_format.txt\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T)
	knitr::kable(exp, \"html\") %>%
  	kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F) %>%
  	kableExtra::footnote(number = c(\"QNAME：比对片段的编号\", 
  	                                \"FLAG：比对情况得分\",
  	                                \"RNAME：比对到基因组上的序列编号\",
  	                                \"POS：比对到基因组上的位置，从1开始计数\",
  	                                \"MAPQ：比对的质量值\",
  	                                \"CIGAR：简要比对信息表达式\",
  	                                \"MRNM：配对片段比对到基因组上的序列编号\",
  	                                \"MPOS：配对片段比对到基因组上的位置\",
  	                                \"ISIZE：插入片段长度\",
  	                                \"SEQ：比对片段序列信息\",
  	                                \"QUAL：比对片段序列的质量信息\",
  	                                \"NH：比对片段比对到基因组上的次数\",
  	                                \"HI：比对片段比对到基因组上的索引\",
  	                                \"AS：局部比对得分\",
  	                                \"nM：成对片段的错配数\")) %>%
  	kableExtra::scroll_box(width = \"800px\", height = \"300px\")
}else{
  print(\"Null Results\")
}
```

<center>表3.1 bam文件格式</center>




### 2.2 比对质控

#### 2.2.1 整体比对情况

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_star_alignment_plot_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.2.3 样本比对情况分布</center>


```{r out.width = \"50%\", echo=FALSE}
table_path = \"$outdir/7.Further_analysis/1.QC/MapQC/Map_stat.txt\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T)
	knitr::kable(exp, \"html\") %>%
  	kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F)
}else{
  print(\"Null Results\")
}
```

<center>表3.1 样本比对情况统计</center>


#### 2.2.2 文库插入片段大小分布

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_qualimap_insert_size_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.2.4 文库插入片段大小分布</center>


#### 2.2.3 GC含量分布

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_qualimap_gc_content_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.2.5 比对上基因组的读长GC含量分布</center>


#### 2.2.4 转录本读长覆盖分布

```{r out.width = \"90%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/Multiqc/multiqc_plots/png/mqc_qualimap_gene_coverage_profile_1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.2.6 转录本读长覆盖分布</center>


## 3. 定量

定量是表达差异显著性分析的基础，目前可在基因水平和转录本水平进行定量。基因水平定量稳健可靠，在算法上容易实现，但无法精细到每个基因的不同转录本。转录本水平定量算法实现起来难度大，精准性不及基因水平定量，但近年来不断有新的算法被开发出来，精准性有显著提高。基因水平定量一般采用简单的reads计数模型，本流程使用RSEM软件对基因和转录本进行定量分析。RSEM分析流程如下图所示：

```{r out.width = \"40%\", echo=FALSE}
img1_path = \"$pipedir/Images/RSEM_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.3.1 RSEM分析流程</center>


### 3.1 定量结果说明

本流程对每个样本分别进行基因水平或转录本水平定量（具体输出格式请参见RSEM软件说明文档），再合并得到所有样本的表达矩阵，第一列为基因或转录本名称，其余列为各样本的表达量（共有两个表，分别为TPM、reads count），如下表所示（TPM矩阵）：

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$outdir/7.Further_analysis/2.Exp_Merge/Exp_TPM.all_0.mat\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T, sep = \"\t\")
	knitr::kable(exp[1:10,], \"html\") %>%
  	kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F) %>%
  	scroll_box(width = \"800px\", height = \"300px\")
}else{
  print(\"Null Results\")
}
```

<center>表3.2 基因表达矩阵（TPM，前10行）</center>


### 3.2基因表达定量质控

#### 3.2.1 基因检测饱和曲线

```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$outdir/7.Further_analysis/1.QC/CountQC/images_GlobalReport/Saturation.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.3.2 基因检测饱和曲线</center>


#### 3.2.2 基因表达分布

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/Exp_TPM.all_0.TPM_dis.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.3.3 基因表达分布图</center>


## 4.差异分析

定量分析完成后，得到所有样本的表达矩阵，便可进行基因或转录本水平的表达差异显著性分析，以便寻找跟实验处理相关的功能基因或转录本。对于有生物学重复实验，本流程采用DESeq2软件进行表达差异显著性分析，将padj小于0.05作为差异显著性标准。对于无生物学重复实验，本流程采用edgeR软件进行表达差异显著性分析，将padj小于0.05，foldchange绝对值大于2作为差异显著性标准。若比较组合中有一组差异基因个数少于100个，则进行调参处理，将pvalue小于0.05作为差异显著性分析标准。表达差异显著性分析包括标准化、离散估计和显著性检验三步。标准化主要是去除测序深度的影响，离散估计能够有效降低假阳性率，而对于无生物学重复实验，无法进行离散估计，需指定一个经验离散值进行处理。实验设计举例如下：


```{r out.width = \"50%\", echo=FALSE}
table_path = \"$pipedir/Images/Normal_design.txt\"
exp = read.table(table_path, header = T)
knitr::kable(exp, \"html\") %>%
  kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F)
```

<center>表3.5 常规实验设计</center>


```{r out.width = \"50%\", echo=FALSE}
table_path = \"$pipedir/Images/Mult_group_design.txt\"
exp = read.table(table_path, header = T)
knitr::kable(exp, \"html\") %>%
  kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F)
```

<center>表3.5 多组样本实验设计</center>


```{r out.width = \"50%\", echo=FALSE}
table_path = \"$pipedir/Images/Pair_design.txt\"
exp = read.table(table_path, header = T)
knitr::kable(exp, \"html\") %>%
  kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F)
```

<center>表3.5 成对样本实验设计</center>


###{.tabset .tabset-dropdown}\n";

open IN, "$outdir/tmp/Design.txt" or die $!;
while(<IN>){
	chomp;
	my $comp_group = $_;
	print "
#### $comp_group

###### 4.1 差异基因列表

根据客户的实验设计，采用相应的统计模型进行表达差异显著性分析，结果如下表所示（分为显著上调基因、显著下调基因以及所有基因三个表格）：

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$outdir/7.Further_analysis/3.DEG/$comp_group.deseq.up_regulate.xls\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T)
	knitr::kable(exp[1:10,], \"html\") %>%
    kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F)
}else{
  print(\"Null Results\")
}
```

<center>表3.3 差异上调基因列表（Top 10）</center>

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$outdir/7.Further_analysis/3.DEG/$comp_group.deseq.down_regulate.xls\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T)
	knitr::kable(exp[1:10,], \"html\") %>%
    	kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F)
}else{
  print(\"Null Results\")
}
```

<center>表3.4 差异基因下调列表（Top 10）</center>


###### 4.2 差异基因注释（目前没有，需要的话后面可以加）

###### 4.3 差异基因火山图

火山图可直观显示表达差异显著性基因的整体分布情况，横坐标表示基因在不同样本中的表达倍数变化(log2FoldChange)，纵坐标表示表达差异的显著性水平(-log10padj)。若比较组合无表达差异显著性基因，默认调整筛选表达差异显著性的阈值进行火山图的绘制。上调基因用红色点表示，下调基因用蓝色点表示，如下图所示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.deseq.Plot-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.4.2 差异基因火山图</center>


###### 4.4 基于差异表达基因的样本聚类

两组以上的实验，可对差异基因集进行聚类分析，将表达模式相近的基因聚在一起，这些基因可能具有共同的功能或参与到共同的代谢途径及信号通路中。本分析流程采用主流的层次聚类方法，对样本的表达值（Normalized reads count）进行聚类，聚类热图如下所示，热图中红色表示高表达，蓝色表示低表达。

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.deseq.Plot-2.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.4.3 差异基因聚类热图</center>


###### 4.5 主成分分析（PCA）


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.deseq.Plot-4.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.4.4 差异基因主成分分析</center>

\n";
}

print "
## 5.功能富集分析

### 5.1通过对差异基因进行富集分析

通过对差异基因进行富集分析，可以找到不同条件下的差异基因与哪些生物学功能或通路显著性相关。本分析流程采用clusterProfiler软件对差异基因集进行GO功能富集分析， KEGG通路富集分析，Reactome通路富集分析。富集分析基于超几何分布原理，如下图所示，其中差异基因集为差异显著分析所得差异基因列表注释到ENTREZ数据库的基因集，背景基因集为所有进行差异显著分析的基因列表注释到ENTREZ数据库的基因集。


```{r out.width = \"50%\", echo=FALSE}
img1_path = \"$pipedir/Images/enrichment_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.1 功能富集分析原理</center>

###{.tabset .tabset-dropdown}\n";

open IN, "$outdir/tmp/Design.txt" or die $!;
while(<IN>){
	chomp;
	my $comp_group = $_;
	print "
#### $comp_group

##### 5.1.1 上调基因功能富集

###### 5.1.1.1 GO功能富集

GO(Gene Ontology)是描述基因功能的综合性数据库，可分为分子功能（Molecular Function），生物过程（biological process）和细胞组成（cellular component）三个部分。本分析流程对显著上、下调基因分别进行富集分析，GO富集以padj小于0.05为显著富集，富集结果如下表所示（分为MF、BP、CC三个表格，不一定都有结果）：



* GO富集结果展示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.GO_BP-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.2 GO生物过程功能富集条形图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.GO_MF-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.2 GO分子功能功能富集条形图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.GO_CC-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.2 GO细胞组分功能富集条形图</center>



```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.GO_BP-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.3 GO生物过程功能富集气泡图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.GO_MF-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.3 GO分子功能功能富集气泡图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.GO_CC-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.3 GO细胞组分功能富集气泡图</center>


###### 5.1.1.2 KEGG通路富集

KEGG(Kyoto Encyclopedia of Genes and Genomes)是整合了基因组、化学和系统功能信息的综合性数据库。本分析流程对显著上、下调基因分别进行富集分析，KEGG富集以padj小于0.05为显著富集，富集结果如下表所示：

* KEGG富集结果展示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.KEGG-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.4 KEGG通路富集条形图</center>


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.KEGG-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.5 KEGG通路富集气泡图举例</center>


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.KEGG.path.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.6 KEGG富集通路图举例</center>


###### 5.1.1.3 Reactome富集分析

Reactome数据库汇集了人类各项反应及生物学通路。本分析流程对显著上、下调基因分别进行富集分析，Reactome富集以padj小于0.05为显著富集，富集结果如下表所示：



* Reactome富集结果展示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.Reactome-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.7 Reactome通路富集条形图</center>


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.up.Reactome-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.8 Reactome通路富集气泡图</center>

##### 5.1.2 下调基因功能富集

###### 5.1.2.1 GO功能富集

GO(Gene Ontology)是描述基因功能的综合性数据库，可分为分子功能（Molecular Function），生物过程（biological process）和细胞组成（cellular component）三个部分。本分析流程对显著上、下调基因分别进行富集分析，GO富集以padj小于0.05为显著富集，富集结果如下表所示（分为MF、BP、CC三个表格，不一定都有结果）：



* GO富集结果展示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.GO_BP-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.9 GO生物过程功能富集条形图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.GO_MF-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.10 GO分子功能功能富集条形图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.GO_CC-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.11 GO细胞组分功能富集条形图</center>



```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.GO_BP-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.12 GO生物过程功能富集气泡图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.GO_MF-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.13 GO分子功能功能富集气泡图</center>

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.GO_CC-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.14 GO细胞组分功能富集气泡图</center>


###### 5.1.2.2 KEGG通路富集

KEGG(Kyoto Encyclopedia of Genes and Genomes)是整合了基因组、化学和系统功能信息的综合性数据库。本分析流程对显著上、下调基因分别进行富集分析，KEGG富集以padj小于0.05为显著富集，富集结果如下表所示：

* KEGG富集结果展示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.KEGG-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.15 KEGG通路富集条形图</center>


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.KEGG-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.16 KEGG通路富集气泡图</center>


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.KEGG.path.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.17 KEGG富集通路图举例</center>


###### 5.1.2.3 Reactome富集分析

Reactome数据库汇集了人类各项反应及生物学通路。本分析流程对显著上、下调基因分别进行富集分析，Reactome富集以padj小于0.05为显著富集，富集结果如下表所示：



* Reactome富集结果展示：

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.Reactome-0.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.18 Reactome通路富集条形图</center>


```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.down.Reactome-1.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.5.19 Reactome通路富集气泡图</center>

\n";
}

print "
### 5.2 GSEA功能富集分析

Gene Set Enrichment Analysis (基因集富集分析)用来评估一个预先定义的基因集的基因在与表型相关度排序的基因表中的分布趋势，从而判断其对表型的贡献。与前面的富集分析不同，GSEA不局限于差异基因，从基因集的富集角度出发，理论上更容易囊括细微但协调性的变化对生物通路的影响。

本流程GSEA功能富集同样包括GO、KEGG以及Reactome三个数据库的功能富集，结果举例如下：

###{.tabset .tabset-dropdown}

\n";

open IN, "$outdir/tmp/Design.txt";
while(<IN>){
	chomp;
	my $comp_group = $_;
	print "
#### $comp_group

##### 5.2.1 GSEA功能富集结果表格

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$outdir/7.Further_analysis/4.GSEA/$comp_group/GSEA.gene_set.report.xls\"
if (file.exists(table_path)){
  gsea_all = read.table(table_path, header = T, sep = \"\t\")
  gsea = gsea_all[gsea_all\$fdr<0.05,]
  if (nrow(gsea)>0){
    knitr::kable(gsea_all[1:10,], \"html\") %>%
      kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F) %>%
      scroll_box(width = \"800px\", height = \"300px\")
  }else{
    print(\"No significant enrichment results\")
  }
}else{
  print(\"No enrichment results\")
}
```

<center>表3.5 GSEA功能富集结果（Top 10）</center>


##### 5.2.2 GSEA功能富集结果可视化

```{r out.width = \"80%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$comp_group.GSEA.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"No significant enrichment results\")
}
```

<center>图3.5.20 GSEA功能富集结果可视化举例</center>

\n";
}

print "
## 6.可变剪接分析

可变剪接是调节基因表达和产生蛋白质组多样性的重要机制，也是RNA-seq的重要分析内容，其形成过程如下图所示：

```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$pipedir/Images/splicing_theory.png\"
if (file.exists(img1_path)){
	knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.6.1 可变剪接形成示意图</center>


本流程使用rMATS软件进行可变剪接分析，主要包括SE、RI、MXE、A5SS、A3SS五种可变剪接事件，如下图所示：

```{r out.width = \"50%\", echo=FALSE}
img1_path = \"$pipedir/Images/splicing_types.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.6.2 rMATS分析可变剪接的五种类型</center>


### 6.1 可变剪接差异分析结果

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$pipedir/Images/Splicing.txt\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T)
	knitr::kable(exp, \"html\") %>%
  	kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F) %>%
  	scroll_box(width = \"800px\", height = \"300px\")
}else{
  print(\"Null Results\")
}
```

<center>表3.6 可变剪接差异分析结果格式举例</center>


### 6.2 差异可变剪接可视化

对于表达差异显著性的可变剪接事件，进行可视化展示，如下图所示。图中跨外显子比对的reads使用连接外显子junction边界的弧线表示。弧线的粗细和比对到junction上的reads数成正比，同时弧线上的数字指出了junction reads的数目。

```{r out.width = \"70%\", echo=FALSE}
img1_path = \"$pipedir/Images/Splicing.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.6.3 差异可变剪接可视化举例</center>


## 7.变异位点分析

变异位点分析是RNA-seq结构分析的重要内容，主要包括先天变异位点和后天体细胞突变位点的检测，对肿瘤等研究具有重要意义。本分析流程使用GATK软件对样本数据进行变异位点分析，并用SnpEff软件对变异位点进行注释，其分析流程如图所示：

```{r out.width = \"60%\", echo=FALSE}
img1_path = \"$pipedir/Images/SNV_detection_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.7.1 GATK分析流程</center>

### 7.1 变异位点检测

变异位点主要分为SNP与INDEL，每个位点的基因型及注释如表7.1所示（只展示了列名及一个位点的注释结果）：

```{r out.width = \"50%\", echo=FALSE}
table_path = \"$pipedir/Images/SNP.txt\"
if (file.exists(table_path)){
	exp = read.table(table_path, header = T, sep = \"\\t\", fill = T)
	knitr::kable(exp, \"html\") %>%
 		kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F) %>%
  	scroll_box(width = \"800px\", height = \"300px\")
}else{
  print(\"Null Results\")
}
```

<center>表3.7 变异位点及注释结果格式举例</center>


## 8.融合基因分析

融合基因是指两个基因的全部或部分序列融合而成的嵌合基因，一般由染色体易位、缺失等原因所致。融合基因首次发现于血液系统的恶性肿瘤中，其中以慢性粒细胞白血病中BCR-ABL的基因融合最为经典，治疗慢性粒细胞白血病的药物伊马替尼/格列卫，其作用靶点就是该融合基因。高通量RNA测序技术因其通量高、成本低、检测精度高和检测范围广等优点大大加快了融合基因的研究，本流程使用STARfusion软件进行融合基因的检测，分析流程如图8.1所示。

```{r out.width = \"50%\", echo=FALSE}
img1_path = \"$pipedir/Images/fusion_detection_theory.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.7.2 STAR-fusion分析流程</center>


融合基因检测结果格式如下表所示：

```{r out.width = \"50%\", echo=FALSE}
library(magrittr)
table_path = \"$pipedir/Images/Fusion_out.txt\"
if (file.exists(table_path)){
  exp = read.table(table_path, header = T)
	knitr::kable(exp, \"html\") %>%
  	kableExtra::kable_styling(bootstrap_options = \"striped\", full_width = F) %>%
  	scroll_box(width = \"800px\", height = \"300px\")
}else{
  print(\"Null Results\")
}
```

<center>表3.8 融合基因检测结果示例</center>


### 各样本融合基因检测结果{.tabset .tabset-dropdown}
\n";

my @sample = `ls $outdir/5.Fusion`;
for my $ss(@sample){
	chomp $ss;
	print"
	
#### $ss

融合基因在染色体上的分布如下图所示，图中的环由多条染色体组成，环内每条线代表一个融合事件，线的两端代表融合事件的断点位置，线的粗细代表融合事件的读长（reads）支持数。

```{r out.width = \"100%\", echo=FALSE}
img1_path = \"$outdir/tmp/Images/$ss.fusion_circos.png\"
if (file.exists(img1_path)){
  knitr::include_graphics(img1_path)
}else{
  print(\"Null Results\")
}
```

<center>图3.7.3 融合基因Circos图</center>

\n";
}

print"

## 9.关键基因临床关联分析

将基因表达水平与临床信息（如肿瘤阶段、病理分型、生存期等），寻找具有生物学意义的关键基因。部分结果举例：

```{r out.width = \"70%\", echo=FALSE}
print(\"Null Results\")
```

\n";
