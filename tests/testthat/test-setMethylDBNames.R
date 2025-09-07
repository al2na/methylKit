context("make sure setting names of known classes works")


# https://github.com/al2na/methylKit/issues/360


data("methylKit")

mdb_diff <- methylDiff.obj |> makeMethylDB(dbdir = tempdir()) 

test_that("check if .setMethylDBNames works", {
  expect_equal(
    mdb_diff |>
      getDBPath() |>
      headTabix() |>
      (function(x) {
        x$V5 <- as.integer(x$V5)
        x$V6 <- as.integer(x$V6)
        return(x)
      })() |>
      .setMethylDBNames() |> names(),
    names(methylDiff.obj)
  )
})
  
mdb_diff |> getDBPath() |> unlink()
