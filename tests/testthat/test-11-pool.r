context("pool checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=methRead( file.list,
            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
            pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(methRead( file.list,
                                  sample.id=list("t1","t2","c1","c2"),assembly="hg18",
                                  pipeline="amp",treatment=c(1,1,0,0),dbtype = "tabix",dbdir="methylDB"))

# unite function
methidh=unite(myobj)
suppressMessages(methidhDB <- unite(mydblist))

test_that("test if output of unite is  methylBase object", {
  expect_is(pool(methidh,sample.ids=c("test","control")), 'methylBase')
  expect_is(pool(methidhDB,sample.ids=c("test","control"),save.db = FALSE), 'methylBase')
})

test_that("test if output of unite is  methylBaseDB object", {
  expect_is(pool(methidhDB,sample.ids=c("test","control")), 'methylBaseDB')
  expect_is(pool(methidh,sample.ids=c("test","control"),save.db = TRUE,dbdir="methylDB"), 'methylBaseDB')
})


test_that("wrong number of sample.ids lead to error", {
  expect_error(pool(methidhDB,sample.ids=c("test","control","onetoomuch")))
  expect_error(pool(methidh,sample.ids=c("test","control","onetoomuch")))
})

unlink("tests/testthat/methylDB",recursive = TRUE)