context("read.bismark checks")


test_that("expect that there is an error when reading an unsorted sam file", {
    expect_that( read.bismark(location=system.file("extdata", "test.fastq_bismark.unsorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18") , 
        throws_error())
})


test_that("check that CpG context can be read-in", {
    expect_that( read.bismark(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18"),
        is_a('methylRaw'))
})

test_that("check that CHG context can be read-in", {
    expect_that( read.bismark(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",read.context="CHG"),
        is_a('methylRaw'))
})

test_that("check that CHH context can be read-in", {
    expect_that( read.bismark(location=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),sample.id="test1",assembly="hg18",read.context="CHH"),
        is_a('methylRaw'))
})
