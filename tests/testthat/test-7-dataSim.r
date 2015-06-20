library("methylKit")
context("dataSim checks")

data("methylKit")

# test dataSim
my.methylBase=dataSim(replicates=4,sample.ids=methylBase.obj@sample.ids,assembly=methylBase.obj@assembly,
                      context=methylBase.obj@context,treatment=methylBase.obj@treatment,
                      destranded=methylBase.obj@destranded,resolution=methylBase.obj@resolution)

test_that("check if dataSim works", {
  expect_that(my.methylBase, 
              is_a('methylBase'))
})
