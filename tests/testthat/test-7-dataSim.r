#library("methylKit")
context("dataSim checks")

data("methylKit")

# test dataSim
my.methylBase=dataSim(replicates=4,sites=20000,treatment=c(1,1,0,0))
my.methylBase2=dataSim(replicates=4,sites=20000,treatment=c(1,1,0,0),effect=c(40,50,60))

test_that("check if dataSim output is a methylBase object", {
  expect_that(my.methylBase, 
              is_a('methylBase'))
})

test_that("check if dataSim output with vectorized effect sizes is a methylBase object", {
  expect_that(my.methylBase2, 
              is_a('methylBase'))
})