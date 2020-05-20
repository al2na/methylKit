context("check if methRead functions work (part 2)")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

dbdir <- "methylDB"


test_that("check if methRead works for amp pipeline", {
  expect_equal(getData(
    methRead(
      file.list[[1]],
      sample.id = "test1",
      assembly = "hg18",
      pipeline = "amp"
    )
  ),
  getData(
    methRead(
      file.list[[1]],
      sample.id = "test1",
      assembly = "hg18",
      pipeline = "amp",
      dbtype = "tabix",
      dbdir = dbdir
    )
  ))
})

test_that("check if methRead works for bismark coverage pipeline", {
  expect_equal(getData(
    methRead(
      location = system.file("extdata", "bismark_coverage_CpG.bed", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      pipeline = "bismarkCoverage",
      mincov = 5
    )
  )
  ,
  getData(
    methRead(
      location = system.file("extdata", "bismark_coverage_CpG.bed", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      pipeline = "bismarkCoverage",
      mincov = 5,
      dbtype = "tabix",
      dbdir = dbdir
    )
  )
  )
})

test_that("check if methRead works for bismark coverage pipeline", {
  expect_equal(getData(
    methRead(
      location = system.file("extdata", "bismark_cytosineReport.txt", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      context = "CpG",
      pipeline = "bismarkCytosineReport",
      mincov = 5
    )
  )
  ,
  getData(
    methRead(
      location = system.file("extdata", "bismark_cytosineReport.txt", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      context = "CpG",
      pipeline = "bismarkCytosineReport",
      mincov = 5,
      dbtype = "tabix",
      dbdir = dbdir
    )
  )
  )
})

test_that("check if methRead works for bismark coverage pipeline", {
  expect_equal(getData(
    methRead(
      location = system.file("extdata", "generic1.CpG.txt",package = "methylKit"),
      pipeline = list(
        fraction = FALSE,
        chr.col = 1,
        start.col = 2,
        end.col = 2,
        coverage.col = 4,
        strand.col = 3,
        freqC.col = 5
      ),
      sample.id = "test1",
      assembly = "hg18"
    ))
  ,
  getData(
    methRead(
      location = system.file("extdata", "generic1.CpG.txt",package = "methylKit"),
      pipeline = list(
        fraction = FALSE,
        chr.col = 1,
        start.col = 2,
        end.col = 2,
        coverage.col = 4,
        strand.col = 3,
        freqC.col = 5
      ),
      sample.id = "test1",
      assembly = "hg18",
      dbtype = "tabix",
      dbdir = dbdir
    )
  )
  )
})


unlink(dbdir,recursive = TRUE)
