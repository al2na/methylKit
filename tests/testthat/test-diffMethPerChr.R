context("check if diffMethPerChr works")

data(methylKit)

dbdir <- "methylDB"

methyldiffDB = makeMethylDB(methylDiff.obj, dbdir = dbdir)

test_that("check if diffMethPerChr works", {
  expect_length(
    diffMethPerChr(
      methylDiff.obj,
      plot = FALSE,
      qvalue.cutoff = 0.01,
      meth.cutoff = 25,
      exclude = NULL
    ),
    2
  )
  expect_equal(names(
    diffMethPerChr(
      methylDiff.obj,
      plot = TRUE,
      qvalue.cutoff = 0.01,
      meth.cutoff = 25,
      exclude = NULL
    )
  )
  , c("rect", "text"))
  expect_warning(
    diffMethPerChr(
      methylDiff.obj,
      plot = TRUE,
      qvalue.cutoff = 0.01,
      meth.cutoff = 25,
      exclude = "chr21"
    )
  )
  expect_length(
    diffMethPerChr(
      methyldiffDB,
      plot = FALSE,
      qvalue.cutoff = 0.01,
      meth.cutoff = 25,
      exclude = NULL
    ),
    2
  )
})
