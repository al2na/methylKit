context("processBismarkAln checks")


test_that("expect that there is an error when reading an unsorted sam file", {
    expect_that( suppressMessages(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.unsorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18")) , 
        throws_error())
})


test_that("check that CpG context can be read-in as methylRaw", {
    expect_is( suppressMessages(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18")),
        'methylRaw')
})

test_that("check that CpG context can be read-in as methylRawDB", {
  expect_is( suppressWarnings(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",save.folder="methylDB",save.db=TRUE)),
     'methylRawDB')
})

test_that("check that CHG context can be read-in as methylraw", {
    expect_is( suppressMessages(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",read.context="CHG")),
        'methylRaw')
})

test_that("check that CHG context can be read-in as methylrawdb", {
  expect_is( suppressWarnings(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",read.context="CHG",save.folder="methylDB",save.db=TRUE)),
             'methylRawDB')
})

test_that("check that CHH context can be read-inas methylraw", {
    expect_is( suppressMessages(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",read.context="CHH")),
        'methylRaw')
})

test_that("check that CHH context can be read-in as methylrawdb", {
  expect_is( suppressWarnings(processBismarkAln(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",read.context="CHH",save.folder="methylDB",save.db=TRUE)),
             'methylRawDB')
})

