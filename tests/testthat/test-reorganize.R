context("check reorganizing")

data("methylKit")

methylBaseDB <- makeMethylDB(methylBase.obj,dbdir = "methylDB")
methylRawListDB_ <- makeMethylDB(methylRawList.obj,dbdir = "methylDB")

## check main reorganize function 

test_that("reorganize works for methylBase", {
  expect_is(reorganize(
    methylBase.obj,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0)
  ),
  "methylBase")
})

test_that("reorganize works for methylBaseDB", {
  expect_is(reorganize(
    methylBaseDB,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0)
  ),
  "methylBaseDB")
})


test_that("reorganize works for methylRawList", {
  expect_is(reorganize(
    methylRawList.obj,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0)
  ),
  "methylRawList")
})

test_that("reorganize works for methylRawListDB", {
  expect_is(reorganize(
    methylRawListDB_,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0)
  ),
  "methylRawListDB")
})

## check reorganize failing
## check if lenghts of sample ids and trtmt don't match

test_that("reorganize fails if lenghts of sample ids and trtmt don't match", {
  expect_error(reorganize(
    methylBase.obj,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0, 0)
  ),
  regexp = "length of sample.ids should be equal to treatment")
})

test_that("reorganize fails if lenghts of sample ids and trtmt don't match", {
  expect_error(reorganize(
    methylBaseDB,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0, 0)
  ),
  regexp = "length of sample.ids should be equal to treatment")
})

test_that("reorganize fails if lenghts of sample ids and trtmt don't match", {
  expect_error(reorganize(
    methylRawList.obj,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0, 0)
  ),
  regexp = "length of sample.ids should be equal to treatment")
})

test_that("reorganize fails if lenghts of sample ids and trtmt don't match", {
  expect_error(reorganize(
    methylRawListDB_,
    sample.ids = c("test1", "ctrl2"),
    treatment = c(1, 0, 0)
  ),
  regexp = "length of sample.ids should be equal to treatment")
})

## check if sample ids not known

test_that("reorganize fails if sample ids not known", {
  expect_error(reorganize(
    methylBase.obj,
    sample.ids = c("test1", "anything"),
    treatment = c(1, 0)
  ),
  regexp = "provided sample.ids is not a subset of the sample ids of the object")
})

test_that("reorganize fails if sample ids not known", {
  expect_error(reorganize(
    methylBaseDB,
    sample.ids = c("test1", "anything"),
    treatment = c(1, 0)
  ),
  regexp = "provided sample.ids is not a subset of the sample ids of the object")
})

test_that("reorganize fails if sample ids not known", {
  expect_error(reorganize(
    methylRawList.obj,
    sample.ids = c("test1", "anything"),
    treatment = c(1, 0)
  ),
  regexp = "provided sample.ids is not a subset of the sample ids of the object")
})

test_that("reorganize fails if sample ids not known", {
  expect_error(reorganize(
    methylRawListDB_,
    sample.ids = c("test1", "anything"),
    treatment = c(1, 0)
  ),
  regexp = "provided sample.ids is not a subset of the sample ids of the object")
})


## check if saves to tabix

test_that("reorganize saves to tabix", {
  expect_is(
    reorganize(
      methylBase.obj,
      sample.ids = c("test1", "ctrl2"),
      treatment = c(1, 0),
      save.db = TRUE,
      dbdir = "methylDB"
    ),
    "methylBaseDB"
  )
})

test_that("reorganize saves to tabix", {
  expect_is(
    reorganize(
      methylRawList.obj,
      sample.ids = c("test1", "ctrl2"),
      treatment = c(1, 0),
      save.db = TRUE,
      dbdir = "methylDB"
    ),
    "methylRawListDB"
  )
})

## check if saves to memory

test_that("reorganize saves to memory", {
  expect_is(
    reorganize(
      methylBaseDB,
      sample.ids = c("test1", "ctrl2"),
      treatment = c(1, 0),
      save.db = FALSE
    ),
    "methylBase"
  )
})

test_that("reorganize saves to memory", {
  expect_is(
    reorganize(
      methylRawListDB_,
      sample.ids = c("test1", "ctrl2"),
      treatment = c(1, 0),
      save.db = FALSE
    ),
    "methylRawList"
  )
})

## check if saves to tabix with suffix

test_that("reorganize saves to tabix with suffix", {
  expect_equal(
    getDBPath(reorganize(
      methylBase.obj,
      sample.ids = c("test1", "ctrl2"),
      treatment = c(1, 0),
      save.db = TRUE,
      dbdir = "methylDB",
      suffix = "new"
    )),
    "methylDB/methylBase_new.txt.bgz")
})

test_that("reorganize saves to tabix with suffix", {
  expect_equal(
    getDBPath(reorganize(
      methylRawList.obj,
      sample.ids = c("test1", "ctrl2"),
      treatment = c(1, 0),
      save.db = TRUE,
      dbdir = "methylDB",
      suffix = "new"
    )),
    c("methylDB/test1_new.txt.bgz","methylDB/ctrl2_new.txt.bgz"))
})

unlink("tests/testthat/methylDB",recursive = TRUE)