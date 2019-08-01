context("unite checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=methRead( file.list,
            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
            pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(methRead( file.list,
                                  sample.id=list("t1","t2","c1","c2"),
                                  assembly="hg18", pipeline="amp",
                                  treatment=c(1,1,0,0), dbtype = "tabix",
                                  dbdir="methylDB"))

test_that("test if output of unite is  methylBase object", {
    expect_is(unite(myobj,destrand=T), 'methylBase')
    expect_is(unite(myobj), 'methylBase')
    expect_is(unite(myobj,min.per.group=1L), 'methylBase')
    expect_is(unite(mydblist,destrand=T,save.db = F), 'methylBase')

})

test_that("test if output of unite is  methylBaseDB object", {
  expect_is(unite(mydblist,destrand=T), 'methylBaseDB')
  expect_is(unite(mydblist), 'methylBaseDB')
  expect_is(unite(mydblist,min.per.group=1L), 'methylBaseDB')
  expect_is(unite(myobj,destrand=T,save.db = T,dbdir="methylDB"), 'methylBaseDB')
})


# test for error when methylRaws are non-overlapping

filtered.obj <- filterByCoverage(methylRawList.obj,
                                 lo.count = 200, 
                                 lo.perc = NULL, 
                                 hi.count = NULL, 
                                 hi.perc = 99.9)
filtered.objDB <- makeMethylDB(filtered.obj, dbdir = "methylDB")

test_that("test if unite stopps if bases overlap.", {
  expect_error(unite(filtered.obj,
                     destrand = TRUE,
                     save.db = TRUE,
                     dbdir = "methylDB",
                     suffix = "filtered_merged"))
  expect_error(unite(filtered.obj,
                     destrand = TRUE))
  expect_error(unite(filtered.objDB,
                     destrand = TRUE,
                     save.db = TRUE,
                     suffix = "filtered_merged"))
  expect_error(unite(filtered.objDB,
                     destrand = TRUE,
                     save.db = FALSE))
})
