context("check if getmethylDiff works")

data(methylKit)

dbdir <- "methylDB"

methyldiffDB = makeMethylDB(methylDiff.obj, dbdir = dbdir)

# check if getmethyldiff works

test_that("test if getMethylDiff works for methylDiff", {
  expect_is(
    getMethylDiff(
      methylDiff.obj,
      difference = 25,
      qvalue = 0.01,
      type = "all"
    ),
    "methylDiff"
  )
  expect_is(
    getMethylDiff(
      methylDiff.obj,
      difference = 25,
      qvalue = 0.01,
      type = "hyper"
    ),
    "methylDiff"
  )
  expect_is(
    getMethylDiff(
      methylDiff.obj,
      difference = 25,
      qvalue = 0.01,
      type = "hypo"
    ),
    "methylDiff"
  )
})

test_that("test if getMethylDiff works for methylDiffDB", {
  expect_is(
    getMethylDiff(
      methyldiffDB,
      difference = 25,
      qvalue = 0.01,
      type = "all"
    ),
    "methylDiffDB"
  )
  expect_is(
    getMethylDiff(
      methyldiffDB,
      difference = 25,
      qvalue = 0.01,
      type = "hyper"
    ),
    "methylDiffDB"
  )
  expect_is(
    getMethylDiff(
      methyldiffDB,
      difference = 25,
      qvalue = 0.01,
      type = "hypo"
    ),
    "methylDiffDB"
  )
})

# check if expected fails work

test_that("test if getMethylDiff fails", {
  expect_error(
    getMethylDiff(
      methylDiff.obj,
      difference = 25,
      qvalue = 0.01,
      type = "any"
    ),
    "Wrong 'type' argument supplied"
  )
  expect_error(
    getMethylDiff(
      methyldiffDB,
      difference = 25,
      qvalue = 0.01,
      type = "any"
    ),
    "Wrong 'type' argument supplied"
  )
})

# check if type conversions work

test_that("test if getMethylDiff saves to tabix", {
  expect_is(
    getMethylDiff(
      methylDiff.obj,
      difference = 25,
      qvalue = 0.01,
      type = "all",
      save.db = TRUE,
      dbdir = dbdir
    ),
    "methylDiffDB"
  )
})

test_that("test if getMethylDiff saves to memory", {
  expect_is(
    getMethylDiff(
      methyldiffDB,
      difference = 25,
      qvalue = 0.01,
      type = "all",
      save.db = FALSE
    ),
    "methylDiff"
  )
})

test_that("test if getMethylDiff saves to tabix with suffix", {
  expect_equal(
    getDBPath(
      getMethylDiff(
      methylDiff.obj,
      difference = 25,
      qvalue = 0.01,
      type = "all",
      save.db = TRUE,
      dbdir = dbdir,
      suffix = "final"
    )),
    "methylDB/methylDiff_final.txt.bgz"
  )
})


unlink(dbdir,recursive = TRUE)