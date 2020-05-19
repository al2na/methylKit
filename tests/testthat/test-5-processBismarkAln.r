context("processBismarkAln checks")

test_that("expect that there is an error when reading an unsorted sam file", {
  expect_that(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.unsorted.min.sam",
                             package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      verbose = FALSE
    ) ,
    throws_error()
  )
  expect_that(
    processBismarkAln(
      location = system.file(
        "extdata",
        "test.fastq_bismark.unsorted_chr.min.sam",
        package = "methylKit"
      ),
      sample.id = "test1",
      assembly = "hg18",
      verbose = FALSE
    ) ,
    throws_error()
  )
})



test_that("check that bismark single-end bam can be read",
          {
            expect_is(
              suppressWarnings(
              processBismarkAln(
                location = system.file("extdata", "test.bismark_single_end.sorted.bam", package = "methylKit"),
                sample.id = "test1",
                assembly = "hg18",
                treatment = NULL,
                verbose = FALSE
              )),
              'methylRaw'
            )
          })

test_that("check that bismark paired-end bam  can be read",
          {
            expect_is(suppressWarnings(
              processBismarkAln(
                location = system.file("extdata", "ctrl.bismark_paired_end.sorted.bam", package = "methylKit"),
                sample.id = "ctrl1",
                assembly = "hg18",
                treatment = NULL,
                verbose = FALSE
            )),
            'methylRaw')
          })

test_that("check that sam file can be read-in", {
  expect_is(suppressMessages(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      verbose = FALSE
    )
  ),
  'methylRaw')
  expect_is(
    suppressMessages(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      save.folder = "methylDB",
      save.db = TRUE,
      verbose = FALSE
  )),
  'methylRawDB')
})

test_that("check that CHG context can be read in", {
  expect_is(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHG",
      verbose = FALSE
  ),
  'methylRaw')
})

test_that("check that CHG context can be read-in as methylrawdb", {
  expect_is(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHG",
      save.folder = "methylDB",
      save.db = TRUE,
      verbose = FALSE
  ),
  'methylRawDB')
})

test_that("check that CHH context can be read-inas methylraw", {
  expect_is(suppressMessages(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHH",
      verbose = FALSE
    )
  ),
  'methylRaw')
})

test_that("check that CHH context can be read-in as methylrawdb", {
  expect_is(suppressWarnings(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHH",
      save.folder = "methylDB",
      save.db = TRUE,
      verbose = FALSE
    )
  ),
  'methylRawDB')
})

# reading multiple files
file.list2=list(system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"),
                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"),
                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"),
                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"))

test_that("check that multiple files can be read-in as methylrawList", {
  expect_is(
    processBismarkAln(
      location = file.list2,
      sample.id = list("test1", "test2", "ctrl1", "ctrl1"),
      assembly = "hg18",
      save.folder = NULL,
      save.context = NULL,
      read.context = "CpG",
      nolap = FALSE,
      mincov = 10,
      minqual = 20,
      phred64 = FALSE,
      treatment = c(1, 1, 0, 0),
      verbose = FALSE
    ),
    'methylRawList'
  )
})

unlink("tests/testthat/methylDB",recursive = TRUE)