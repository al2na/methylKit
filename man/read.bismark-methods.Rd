\docType{methods}
\name{read.bismark}
\alias{read.bismark}
\alias{read.bismark,character,character,character-method}
\alias{read.bismark,list,list,character-method}
\title{Read from sorted Bismark SAM files}
\usage{
  read.bismark(location, sample.id, assembly,
    save.folder = NULL, save.context = c("CpG"),
    read.context = "CpG", nolap = FALSE, mincov = 10,
    minqual = 20, phred64 = FALSE, treatment)
}
\arguments{
  \item{location}{location of sam file(s). If multiple
  files are given this argument must be a list.}

  \item{sample.id}{the id(s) of samples in the same order
  as file.  If multiple sam files are given this arugment
  must be a list.}

  \item{save.folder}{The folder which will be used to save
  methylation call files, if set to NULL no methylation
  call file will be saved as a text file.  The files saved
  can be read into R in less time using \code{read}
  function in \code{methylKit}}

  \item{save.context}{A character vector consisting
  following strings: "CpG","CHG","CHH".  The methylation
  percentages for these methylation contexts will be saved
  to save.folder}

  \item{read.context}{One of the 'CpG','CHG','CHH' or
  'none' strings.  Determines what type of methylation
  context will be read-in to the memory which can be
  immediately used for analysis.  If given as 'none',
  read.bismark will not return any object, but if a
  save.folder argument given it will save the methylation
  percentage call files.}

  \item{assembly}{string that determines the genome
  assembly. Ex: mm9,hg18 etc.}

  \item{nolap}{if set to TRUE and the SAM file has
  paired-end reads, the one read of the overlapping
  paired-end read pair will be ignored for methylation
  calling.}

  \item{mincov}{minimum read coverage to call a methylation
  status for a base.}

  \item{minqual}{minimum phred quality score to call a
  methylation status for a base.}

  \item{phred64}{logical ( default: FALSE) you will not
  need to set this TRUE, Currently bismark gives only
  phred33 scale}

  \item{treatment}{treatment vector only to be used when
  location and sample.id parameters are \code{list}s and
  you are trying to read-in multiple samples that are
  related to eachother in down-stream analysis.}
}
\value{
  \code{methylRaw} or \code{methylRawList} object
}
\description{
  The function calls methylation percentage per base from
  sorted Bismark SAM files and reads methylation
  information as \code{methylRaw} or \code{methylRawList}
  object. Bismark is a popular aligner for high-throughput
  bisulfite sequencing experiments and it outputs its
  results in SAM format by default. Bismark SAM format
  contains aligner specific tags which are absolutely
  necessary for methylation percentage calling. SAM files
  from other aligners will not work with this function.
}
\note{
  SAM files should be sorted with samtools sort or unix
  sort. Other sorting methods can alter the order of
  fields(columns) in the SAM file and that will result in
  an error when using \code{read.bismark()}.
}
\examples{
# reading one bismark file:
my.file=system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                                                       package = "methylKit")
obj=read.bismark(my.file,"test",assembly="hg18",save.folder=NULL,
                 save.context="CpG",read.context="CpG")

# reading multiple files
file.list2=list(system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
               system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
              system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
               system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"))

 objs=read.bismark(location=file.list2
             ,sample.id=list("test1","test2","ctrl1","ctrl1"),assembly="hg18",
             save.folder=NULL,save.context=NULL,read.context="CpG",
             nolap=FALSE,mincov=10,minqual=20,phred64=FALSE,treatment=c(1,1,0,0))
}

