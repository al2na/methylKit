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



## dummy data containing NA in the counts
missingData <- 
  "chrBase            chr       base     strand   coverage          freqC    freqT
chr1.662657      chr1     662657 F          53        100       0
chr1.662678      chr1     662678 F          106       1.88679245283019        98.1132075471698
chr1.662692      chr1     662692 F          106       100       0
chr1.662694      chr1     662694 F          106       100       0
chr1.680125      chr1     680125 F          NA       NA       NA
chr1.680466      chr1     680466 F          NA       NA       NA
chr1.714259      chr1     714259 R          NA       NA       NA
chr1.714262      chr1     714262 R          NA       NA       NA
chr1.714265      chr1     714265 R          NA       NA       NA
chr1.714278      chr1     714278 R          NA       NA       NA
chr1.714293      chr1     714293 R          NA       NA       NA
chr1.714544      chr1     714544 R          178       0          100
chr1.714547      chr1     714547 R          178       78.6516853932584        21.3483146067416
chr1.714566      chr1     714566 R          178       0.561797752808989      99.438202247191"

## dump dummy data in file
write(missingData,file = "missingData.txt")
## clean up
on.exit(unlink("missingData.txt"))
issue <- "missingData.txt"

# read the files to a methylRawList object: myobj
## NA's are propagated, which should not happen
test_that("check if methRead stops with NA in Input", {
  expect_warning(methRead(issue,sample.id="test",assembly="hg19",
                         context="CpG",header = TRUE,sep = " ",
                         mincov = 0))
})
# --> leads to same issue 
test_that("check if methRead to tabix stops with NA in Input", {
  expect_warning(methRead(issue,sample.id="test",assembly="hg19",
                        context="CpG",header = TRUE,sep = " ",
                        mincov = 0,dbtype = "tabix"))
})


