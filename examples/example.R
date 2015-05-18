#source("~/Dropbox/PAPERS/R-devel/methylKit/R/backbone.R")
#source("~/Dropbox/PAPERS/R-devel/methylKit/R/diffMeth.R")
#source("~/Dropbox/PAPERS/R-devel/methylKit/R/annotate.R")

# R CMD BUILD --no-resave-data --keep-empty-dirs --install-args=install-tests methylKit
# R CMD INSTALL --install-tests --example methylKit_0.1.tar.gz

#file.list=list("~/Dropbox/PAPERS/R-devel/methylKit/data/test1.myCpG.txt",
#                "~/Dropbox/PAPERS/R-devel/methylKit/data/test2.myCpG.txt",
#                "~/Dropbox/PAPERS/R-devel/methylKit/data/control1.myCpG.txt",
#                "~/Dropbox/PAPERS/R-devel/methylKit/data/control2.myCpG.txt")

library(methylKit)


file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list, sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

#-----------------------------------------------------------------------------------------------

# read 5hmC and 5mC files
hmc.file=system.file("extdata", "test1.myCpG.txt", package = "methylKit")
mc.file =system.file("extdata", "test2.myCpG.txt", package = "methylKit")
my5hmC=read( hmc.file,sample.id="hmc",assembly="hg18")
my5mC =modRead( mc.file,sample.id="mc",assembly="hg18")
# adjusting the 5mC levels using 5hmC levels
adjusted.5mC=adjust.methylC(my5mC,my5hmC)
adjusted.6mC=adjustMethylC(my5mC,my5hmC)

my.file=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit")
obj=read.bismark(my.file,"test",assembly="hg18",save.folder=NULL, save.context="CpG",read.context="CpG")

data(methylKit)
cpg.gr=read.bed(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),remove.unsual=TRUE)
annotate.WithFeature(methylDiff.obj,cpg.gr,strand=FALSE,extend=0,feature.name="CpGi")

cpg.obj=read.feature.flank(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),feature.flank.name=c("CpGi","shores"))
annotate.WithFeature.Flank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="Shores")

gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
annotate.WithGenicParts(methylDiff.obj, gene.obj)

# get differentially methylated bases/regions with specific cutoffs
all.diff=get.methylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="all")
# get hyper-methylated
hyper=get.methylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="hyper")
# get hypo-methylated
hypo=getMethylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="hypo")

my.loc=system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit")
cpg.obj=read.feature.flank(location=my.loc,feature.flank.name=c("CpGi","shores"))

gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))


#--------------------------------------------------------------------------------------------------------------


# get methylation statistics on second file "test2" in myobj which is a class of methylRawList
getMethylationStats(myobj[[2]],plot=F,both.strands=F)

# plot methylation statistics on second file "test2" in myobj which is a class of methylRawList
getMethylationStats(myobj[[2]],plot=T,both.strands=F)

# see what the data looks like for sample 2 in myobj methylRawList
head(myobj[[2]])

# merge all samples to one table by using base-pair locations that are covered in all samples
# setting destrand=T, will merge reads on both strans of a CpG dinucleotide. This provides better 
# coverage, but only advised when looking at CpG methylation
meth=unite(myobj,destrand=T)

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

# calculate differential methylation
myDiff=calculateDiffMeth(meth)

# check how data part of methylDiff object looks like
head( myDiff )


# read-in transcript locations to be used in annotation
gene.obj=read.transcript.features(system.file("tests", "refseq.hg18.bed.txt", package = "methylKit"))

# annotate differentially methylated Cs with promoter/exon/intron using annotation data
annotate.WithGenicParts(myDiff,gene.obj)
