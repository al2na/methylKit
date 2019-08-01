context("test adjustMethylC check")


# read 5hmC and 5mC files
hmc.file=system.file("extdata", "test1.myCpG.txt", package = "methylKit")
mc.file =system.file("extdata", "test2.myCpG.txt", package = "methylKit")

my5hmC=methRead( hmc.file,sample.id="hmc",assembly="hg18")
my5mC =methRead( mc.file,sample.id="mc",assembly="hg18")
my5hmCList = methylRawList( my5hmC,my5hmC, treatment = c(1,1))
my5mCList = methylRawList( my5mC,my5mC, treatment = c(1,1))

# adjusting the 5mC levels using 5hmC levels

test_that("test if output of adjustMethylC is  methylRaw or methylRawList object", {
  expect_is(adjustMethylC(my5mC,my5hmC), 'methylRaw')
  expect_is(adjustMethylC(my5mCList,my5hmCList), 'methylRawList')
  
})

makeMethylDB(my5mCList,dbdir = "methylDB")


suppressMessages(expr = {
  my5mCListDB = makeMethylDB(my5mCList)
  my5hmCListDB = makeMethylDB(my5hmCList)
  })
test_that("test if output of adjustMethylC is  methylRawDB object", {
  expect_is(adjustMethylC(my5mC,my5hmC,save.db = TRUE,dbtype="tabix",dbdir="methylDB"), 'methylRawDB')
  expect_is(adjustMethylC(makeMethylDB(my5mC),makeMethylDB(my5hmC)), 'methylRawDB')
  expect_is(adjustMethylC(my5mCListDB,my5hmCListDB), 'methylRawListDB')
})


unlink("tests/testthat/methylDB*",recursive = TRUE)