

<a name="logo"/>
<div align="center">
<img src="https://dl.dropboxusercontent.com/u/1373164/methylKit_logo.png" alt="methylKit Logo"  ></img>
</a>
</div>

methylKit 
========

Build Status [![Build Status](https://travis-ci.org/al2na/methylKit.svg?branch=development)](https://travis-ci.org/al2na/methylKit)
[![GitHub release](https://img.shields.io/github/release/al2na/methylKit.svg)](https://github.com/al2na/methylKit/releases)


# Introduction

*methylKit* is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) package 
for DNA methylation analysis and annotation from high-throughput bisulfite sequencing. The 
package is designed to deal with sequencing data from 
[RRBS](http://www.nature.com/nprot/journal/v6/n4/abs/nprot.2010.190.html) and its variants,
but also target-capture methods such as [Agilent SureSelect 
methyl-seq](http://www.halogenomics.com/sureselect/methyl-seq). 
In addition, methylKit can 
deal with base-pair resolution data for 5hmC obtained from Tab-seq or oxBS-seq. It can also 
handle whole-genome bisulfite sequencing data if proper input format is provided.

## Current Features

 * Coverage statistics
 * Methylation statistics
 * Sample correlation and clustering
 * Differential methylation analysis 
 * Feature annotation and accessor/coercion functions 
 * Multiple visualization options  
 * Regional and tiling windows analysis
 * (Almost) proper [documentation](http://rpubs.com/al2na/methylKit)
 * Reading methylation calls directly from [Bismark(Bowtie/Bowtie2](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/) alignment files
 * Batch effect control
 * Multithreading support (for faster differential methylation calculations) 
 * Coercion to objects from Bioconductor package GenomicRanges
 * Reading methylation percentage data from generic text files





## Staying up-to-date

You can subscribe to our googlegroups page to get the latest information about new releases and features (low-frequency, only updates are posted)

- https://groups.google.com/forum/#!forum/methylkit

To ask questions please use methylKit_discussion forum

- https://groups.google.com/forum/#!forum/methylkit_discussion

You can also check out the blogposts we make on using methylKit

- http://zvfak.blogspot.de/search/label/methylKit

-------

## Installation

in R console,
```r
library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
  repos=BiocInstaller::biocinstallRepos(),
  dependencies=TRUE)
```
if this doesn't work, you might need to add `type="source"` argument.

### Install the development version
```r
library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
  repos=BiocInstaller::biocinstallRepos(),ref="development",
  dependencies=TRUE)
```
if this doesn't work, you might need to add `type="source"` argument.


-------

# How to Use

Typically, bisulfite converted reads are aligned to the genome and % methylation value per base is calculated by processing alignments. *`methylKit`* takes that  % methylation value per base information as input. Such input file may be obtained from [AMP pipeline](http://code.google.com/p/amp-errbs/) for aligning RRBS reads. A typical input file looks like this:

```
chrBase	chr	base	strand	coverage	freqC	freqT
chr21.9764539	chr21	9764539	R	12	25.00	75.00
chr21.9764513	chr21	9764513	R	12	0.00	100.00
chr21.9820622	chr21	9820622	F	13	0.00	100.00
chr21.9837545	chr21	9837545	F	11	0.00	100.00
chr21.9849022	chr21	9849022	F	124	72.58	27.42
chr21.9853326	chr21	9853326	F	17	70.59	29.41

```


*`methylKit`* reads in those files and performs basic statistical analysis and annotation for differentially methylated regions/bases. Also a tab separated text file with a generic format can be read in, such as methylation ratio files from [BSMAP](http://code.google.com/p/bsmap/), see [here](http://zvfak.blogspot.com/2012/10/how-to-read-bsmap-methylation-ratio.html) for an example. Alternatively, `read.bismark` function can read SAM file(s) output by [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/)(using bowtie/bowtie2) aligner (the SAM file must be sorted based on chromosome and read start). The sorting must be done by unix sort or samtools, sorting using other tools may change the column order of the SAM file and that will cause an error. 

Below, there are several options showing how to do basic analysis with *`methylKit`*.

## Documentation ##
 * You can look at the vignette [here](http://rpubs.com/al2na/methylKit)
 * You can check out the [slides](http://methylkit.googlecode.com/files/methylKitTutorialSlides_2013.pdf ) for a tutorial at EpiWorkshop 2013
 * You can check out the [tutorial](http://methylkit.googlecode.com/files/methylKitTutorial_feb2012.pdf) prepared for  EpiWorkshop 2012
 * You can see the code snippet below


## Example analysis ##

### Read methylation files ###
```r
library(methylKit)


# this is a list of example files, ships with the package
# for your own analysis you will just need to provide set of paths to files
# you will not need the "system.file(..."  part
file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawList object: myobj
myobj=read( file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",treatment=c(1,1,0,0))

```

### Get descriptive stats on methylation ###
```r
# get methylation statistics on second file "test2" in myobj which is a class of methylRawList
getMethylationStats(myobj[[2]],plot=F,both.strands=F)

# plot methylation statistics on second file "test2" in myobj which is a class of methylRawList
getMethylationStats(myobj[[2]],plot=T,both.strands=F)
```

### Get bases covered by all samples and cluster samples ###
```r
# see what the data looks like for sample 2 in myobj methylRawList
head(myobj[[2]])

# merge all samples to one table by using base-pair locations that are covered in all samples
# setting destrand=TRUE, will merge reads on both strans of a CpG dinucleotide. This provides better 
# coverage, but only advised when looking at CpG methylation
meth=unite(myobj,destrand=TRUE)

# merge all samples to one table by using base-pair locations that are covered in all samples
# 
meth=unite(myobj)

# cluster all samples using correlation distance and return a tree object for plclust
hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

# cluster all samples using correlation distance and plot hiarachical clustering
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

# screeplot of principal component analysis.
PCASamples(meth, screeplot=TRUE)

# principal component anlaysis of all samples.
PCASamples(meth)
```
### Calculate differential methylation ###
Before differential methylation calculation, consider filtering high coverage bases to remove potential PCR bias using `filterByCoverage()`. In addition, consider normalizing read coverages between samples to avoid bias introduced by systematically more sequenced samples, using `normalizeCoverage()`.

```r
# calculate differential methylation p-values and q-values
myDiff=calculateDiffMeth(meth)

# check how data part of methylDiff object looks like
head( myDiff )

# get differentially methylated regions with 25% difference and qvalue<0.01
myDiff25p=get.methylDiff(myDiff,difference=25,qvalue=0.01)

# get differentially hypo methylated regions with 25% difference and qvalue<0.01
myDiff25pHypo =get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hypo") 

# get differentially hyper methylated regions with 25% difference and qvalue<0.01
myDiff25pHyper=get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
```

### Annotate differentially methylated bases/regions ###
```r
# read-in transcript locations to be used in annotation
# IMPORTANT: annotation files that come with the package (version >=0.5) are a subset of full annotation
# files. Download appropriate annotation files from UCSC (or other sources) in BED format
library(genomation) # install from BioC
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))

# annotate differentially methylated Cs with promoter/exon/intron using annotation data
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
```

SEE PACKAGE VIGNETTE and TUTORIAL (both hyper-linked above) FOR MORE

## Downloading Annotation Files
You can download annotation files from UCSC table browser for your genome of interest. Go to  [http://genome.ucsc.edu/cgi-bin/hgGateway]. On the top menu click on "tools" then "table browser". Select your "genome" of interest and "assembly" of interest from the drop down menus. Make sure you select the correct genome and assembly. Selecting wrong genome and/or assembly will return unintelligible results in downstream analysis. 

From here on you can either download *gene annotation* or *CpG island annotation*.

1. For gene annotation, select _"Genes and Gene prediction tracks"_ from the *"group"* drop-down menu. Following that, select _"Refseq Genes"_ from the *"track"* drop-down menu. Select _"BED- browser extensible data"_ for the *"output format"*. Click *"get output"* and on the following page click *"get BED"* without changing any options. save the output as a text file.
2. For CpG island annotation, select _"Regulation"_ from the *"group"* drop-down menu. Following that, select _"CpG islands"_ from the *"track"* drop-down menu. Select _"BED- browser extensible data"_  for the *"output format"*. Click *"get output"* and on the following page click *"get BED"* without changing any options. save the output as a text file.


In addition, you can check this tutorial to learn how to download any track from UCSC in BED format (http://www.openhelix.com/cgi/tutorialInfo.cgi?id=28)




-------
# R script for Genome Biology publication
The most recent version of the R script in the Genome Biology manuscript is [here](http://code.google.com/p/methylkit/downloads/list?q=label:AdditionalFile4 ).

-------
# Citing methylKit
If you used methylKit please cite:


 * Altuna Akalin, Matthias Kormaksson, Sheng Li, Francine E. Garrett-Bakelman, Maria E. Figueroa, Ari Melnick, Christopher E. Mason. _(2012)_. *"[methylKit: A comprehensive R package for the analysis of genome-wide DNA methylation profiles.](http://genomebiology.com/2012/13/10/R87/)"* _Genome Biology_ , 13:R87.


and also consider citing the following publication as a use-case with specific cutoffs:

 * Altuna Akalin, Francine E. Garrett-Bakelman, Matthias Kormaksson, Jennifer Busuttil, Lu Zhang, Irina Khrebtukova, Thomas A. Milne, Yongsheng Huang, Debabrata Biswas, Jay L. Hess, C. David Allis, Robert G. Roeder, Peter J. M. Valk, Bob Löwenberg, Ruud Delwel, Hugo F. Fernandez, Elisabeth Paietta, Martin S. Tallman, Gary P. Schroth, Christopher E. Mason, Ari Melnick, Maria E. Figueroa. _(2012)_. *"[Base-Pair Resolution DNA Methylation Sequencing Reveals Profoundly Divergent Epigenetic Landscapes in Acute Myeloid Leukemia.](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002781)"* _PLoS Genetics_ 8(6).

-------
# Contact & Questions
e-mail to [methylkit_discussion@googlegroups.com](mailto:methylkit_discussion@googlegroups.com ) or post a question using [the web interface](https://groups.google.com/forum/#!forum/methylkit_discussion).

if you are going to submit bug reports or ask questions, please send sessionInfo() output from R console as well.

Questions are very welcome, although we suggest you read the paper, documentation(function help pages and the vignette) and [ blog entries](http://zvfak.blogspot.com/search/label/methylKit) first. The answer to your question might be there already.

-------
# Contribute to the development
See the [trello board](https://trello.com/b/k2kv1Od7/methylkit) for methylKit development. You can contribute to the methylKit development via github ([http://github.com/al2na/methylKit/]) by checking out "development" branch, making your changes and doing a pull request (all of these should be done on the "development" branch NOT  on the "master" branch). In addition, you should:

 * Bump up the version in the DESCRIPTION file on the 4th number. For example, the master branch has the version numbering as in "X.Y.Z". If you make a change to the development branch you should bump up the version in the DESCRIPTION file to "X.Y.Z.1". If there are already changes done in the development just bump up the fourth number.

 * Add your changes to the NEWS file as well under the correct version and appropriate section. Attribute the changes to yourself, such as "Contributed by X"
 
License
---------
Artistic License/GPL
