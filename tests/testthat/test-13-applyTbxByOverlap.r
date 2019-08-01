context("applyTbxByOverlap checks")
data(methylKit)

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

methylRawListDB.obj=methRead(file.list,
                             sample.id=list("test1","test2","ctrl1","ctrl2"),
                             assembly="hg18",treatment=c(1,1,0,0),
                             dbtype = "tabix",dbdir = "methylDB")

methylBaseDB.obj=unite(methylRawListDB.obj)

methylDiffDB.obj = calculateDiffMeth(methylBaseDB.obj)

# define the windows of interest as a GRanges object, this can be any set
# of genomic locations
library(GenomicRanges)
my.win=GRanges(seqnames="chr21",
               ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )

res.df <- getTabixByOverlap(tbxFile = methylBaseDB.obj@dbpath,granges = my.win,return.type = "data.frame")
res.dt <- getTabixByOverlap(tbxFile = methylBaseDB.obj@dbpath,granges = my.win,return.type = "data.table")

f <- function(x) return(x)

tmp1=applyTbxByOverlap(tbxFile = methylBaseDB.obj@dbpath,
                      ranges=my.win,
                      chunk.size=10,
                      return.type="data.table",FUN = f) 

tmp2=applyTbxByOverlap(tbxFile = methylBaseDB.obj@dbpath,
                      ranges=my.win,
                      chunk.size=10,
                      return.type="data.frame",FUN = f) 

tmp3=applyTbxByOverlap(tbxFile = methylBaseDB.obj@dbpath,
                      ranges=my.win,
                      chunk.size=10,
                      dir = dirname(methylBaseDB.obj@dbpath),
                      filename = paste0(basename(methylBaseDB.obj@dbpath),".out"),
                      return.type="tabix",FUN = f) 


test_that("check if applyTbxByOverlap returns expected results for different return types." ,{
  expect_identical(tmp1,res.dt)
  expect_identical(tmp2,res.df)
  expect_identical(headTabix(tmp3,nrow = nrow(res.dt),return.type = "data.table"),res.dt)
})

unlink("tests/testthat/methylDB",recursive = TRUE)