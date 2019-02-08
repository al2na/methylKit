context("test file list and methRead and getMethylationStats check")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=methRead( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),
		assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(
	methRead( file.list,sample.id=list("test1","test2","ctrl1","ctrl2"),
	assembly="hg18",pipeline="amp",treatment=c(1,1,0,0),
	dbtype = "tabix",dbdir="methylDB"))

mydblist2 = suppressMessages(
	methRead( file.list,sample.id=list("test1","test2","ctrl1","ctrl2"),
	assembly="hg18",pipeline="amp",treatment=c(1,1,0,0),
	dbtype = "tabix"))


mydb = suppressMessages(
	methRead( mydblist[[1]]@dbpath,sample.id="test1",
	assembly="hg18",dbtype = "tabix",dbdir="methylDB"))


test_that("check if there are 4 test files in the file.list",{
    expect_equal(length(file.list),4)
})

test_that("if methRead return a methylRawlist", {
    expect_is(myobj, 'methylRawList')
})

test_that("if methRead return a methylRawListDB", {
  expect_is(mydblist, 'methylRawListDB')
})

test_that("if methRead without given dir return same methylRawListDB", {
  expect_identical(as(mydblist,"methylRawList"),as(mydblist2,"methylRawList"))
})

test_that("if methRead of database return a methylRawDB", {
  expect_is(mydb, 'methylRawDB')
})

test_that("getMethylationStats on methylRawDB works", {
  expect_output(getMethylationStats(mydblist[[2]],plot=F,both.strands=F),
                'methylation statistics per base')
})

test_that("getMethylationStats on methylRaw works", {
    expect_output(getMethylationStats(myobj[[2]],plot=F,both.strands=F),
                  'methylation statistics per base')
})


