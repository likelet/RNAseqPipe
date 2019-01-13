# SYSUCC-RNAseqPipe installation 
> adjusted from nf-core  

- [Installation](#installation)
  * [1) Install NextFlow](#1--install-nextFlow)
  * [2) Install the pipeline](#2--install-the-pipeline)
  * [3) Pipeline configuration](#3--pipeline-configuration)
      - [3.1) Software deps: Docker](#31--software-deps--docker)
      - [3.2) Software deps: Singularity](#32--software-deps--singularity)
      - [3.3) Software deps: conda](#33--software-deps--conda) 
  * [4) Install Software Manually](#4--install-software-manually)


## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

See [nextflow.io](https://www.nextflow.io/) for further instructions on how to install and configure Nextflow.

## 2) Install the pipeline  

#### 2.1) Development  

The whole pipeline is hosted in github. To install, you'll need to download and transfer the pipeline files manually:

```bash
mkdir -p ~/my-pipelines/
git clone https://github.com/likelet/cirPipe.git 
cd ~/my_data/
nextflow run ~/my-pipelines/circPipe
```

To stop nextflow from looking for updates online, you can tell it to run in offline mode by specifying the following environment variable in your ~/.bashrc file:

```bash
export NXF_OFFLINE='TRUE'
```

#### 2.2) Development

If you would like to make changes to the pipeline, it's best to make a fork on GitHub and then clone the files. Once cloned you can run the pipeline directly as above.


## 3) Pipeline configuration
By default, the pipeline runs with the `standard` configuration profile. This uses a number of sensible defaults for process requirements and is suitable for running on a simple (if powerful!) basic server. You can see this configuration in [`conf/base.config`](../conf/base.config).

Be warned of two important points about this default configuration:

1. The default profile uses the `local` executor
    * All jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.
    * See the [nextflow docs](https://www.nextflow.io/docs/latest/executor.html) for information about running with other hardware backends. Most job scheduler systems are natively supported.
2. Nextflow will expect all software to be installed and available on the `PATH`

#### 3.1) Software deps: Docker
First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, running the pipeline with the option `-profile standard,docker` tells Nextflow to enable Docker for this run. An image containing all of the software requirements will be automatically fetched and used from dockerhub (https://hub.docker.com/r/likelet/cirpipe/).

#### 3.2) Software deps: Singularity
If you're not able to use Docker then [Singularity](http://singularity.lbl.gov/) is a great alternative.
The process is very similar: running the pipeline with the option `-profile standard,singularity` tells Nextflow to enable singularity for this run. An image containing all of the software requirements will be automatically fetched and used from singularity hub.

If running offline with Singularity, you'll need to download and transfer the Singularity image first:

```bash
singularity pull --name circPipe.simg shub://likelet/cirpipe
```

Once transferred, use `-with-singularity` and specify the path to the image file:

```bash
nextflow run /path/to/circPipe -with-singularity circPipe.simg
```

Remember to pull updated versions of the singularity image if you update the pipeline.


#### 3.3) Software deps: conda
If you're not able to use Docker _or_ Singularity, you can instead use conda to manage the software requirements.
This is slower and less reproducible than the above, but is still better than having to install all requirements yourself!
The pipeline ships with a conda environment file and nextflow has built-in support for this.
To use it first ensure that you have conda installed (we recommend [miniconda](https://conda.io/miniconda.html)), then follow the same pattern as above and use the flag `-profile standard,conda`

## 4) Install Software Manually
If you want to install all the softwares manually in local, you can follow the steps:
First you should ensure that you have conda installed (we recommend [miniconda](https://conda.io/miniconda.html)), then you can use the `conda install` to install the necessary softwares.
```bash
#install the softwares and dependencies
conda install bwa
conda install circexplorer2
conda install star
conda install multiqc
conda install samtools
conda install bowtie
conda install bowtie2
conda install fastp
conda install perl
```
For [find_circ](https://github.com/marvin-jens/find_circ), you should create an environment based on python2.7 and install some dependencies.
```bash
#install the dependencies
conda create -n tools_in_python2 python=2.7
source activate tools_in_python2
conda install pysam
conda install numpy
source deactivate
#install the find_circ
git clone http://github.com/marvin-jens/find_circ.git
```
For [KNIFE](https://github.com/lindaszabo/KNIFE), you should create an environment based on python2.7 and install some dependencies.
```bash
#install the dependencies
source activate tools_in_python2
conda install bowtie
conda install bowtie2
conda install samtools
conda install biopython
conda install bcbiogff
source deactivate
#install the KNIFE
wget https://github.com/lindaszabo/KNIFE/archive/v1.4.tar.gz
tar zxvf v1.4.tar.gz
rm v1.4.tar.gz
mv KNIFE-1.4 KNIFE
cd KNIFE/circularRNApipeline_Standalone/analysis
chmod a+x *
```
For [Mapsplice](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), you should create an environment based on python3.6 and install some dependencies.
```bash
#install the mapsplice and dependencies
conda create -n tools_in_python3 python=3.6
source activate tools_in_python3
conda install bowtie
conda install mapsplice
source deactivate
```
For [Segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/), you should create an environment based on python3.6 and install some dependencies.
```bash
#install the segemehl and dependencies
source activate tools_in_python3
conda install segemehl
source deactivate
```
For [CIRI](http://sourceforge.net/projects/ciri), you should follow the commands:
```bash
#install the ciri
wget http://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip
unzip CIRI-full_v2.0.zip
rm CIRI-full_v2.0.zip
mv CIRI-full_v2.0 CIRI
```
For the analysis processes using R scripts, you should install the following packages:
```bash
#install the R packages
conda install r-ggplot2
conda install r-pheatmap
conda install bioconductor-edger
conda install bioconductor-deseq2
conda install bioconductor-limma
conda install bioconductor-chipseeker
conda install bioconductor-fgsea
conda install bioconductor-bsgenome.hsapiens.ucsc.hg38
conda install bioconductor-clusterprofiler
conda install bioconductor-reactomepa
conda install bioconductor-pathview
conda install r-dplyr
conda install r-readr
conda install r-rmarkdown
conda install r-ggsci
conda install r-circlize
conda install r-rcolorbrewer
conda install r-venndiagram
```
#### Build necesssary Index
As to build the different indexes, you can simply download from [igenome](https://support.illumina.com/sequencing/sequencing_software/igenome.html), or using your own genome file and annotation file with the command lines introduced in `Parameters` section.

## Appendices

