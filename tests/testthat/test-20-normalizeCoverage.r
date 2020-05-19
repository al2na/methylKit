context("test normalizeCoverage")

data("methylKit")


dbdir <- "methylDB"
mydblist <- makeMethylDB(methylRawList.obj,dbdir = dbdir)


test_that("normalizeCoverage works", {
  # normalize by the median coverage
  expect_is(normalizeCoverage(methylRawList.obj, method = "median"),
            'methylRawList')
  # normalize by mean coverage and save to database in folder methylDB
  expect_is(
    normalizeCoverage(
      methylRawList.obj,
      method = "mean",
      save.db = TRUE,
      dbdir = dbdir
    ),
    'methylRawListDB'
  )
  # normalize by the median coverage
  expect_is(normalizeCoverage(mydblist, method = "median"),
            'methylRawListDB')
  # normalize by mean coverage and save to database in folder methylDB
  expect_is(
    normalizeCoverage(
      mydblist,
      method = "mean",
      save.db = FALSE,
      dbdir = dbdir
    ),
    'methylRawList'
  )
  
  expect_error(
    normalizeCoverage(
      methylRawList.obj,
      method = "max"
    )
  )
  expect_error(
    normalizeCoverage(
      mydblist,
      method = "max"
    )
  )
  
  
})


unlink(dbdir,recursive = TRUE)
