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

data(methylKit)

filtered.obj <- filterByCoverage(methylRawList.obj,
                                 lo.count = 200, 
                                 lo.perc = NULL, 
                                 hi.count = NULL, 
                                 hi.perc = 99.9)
filtered.objDB <- makeMethylDB(filtered.obj, dbdir = "methylDB")
on.exit(unlink("methylDB",recursive = TRUE), add = TRUE, after = TRUE)

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


test_that("test unite with destranding and different resolution and strands", {
  data(methylKit)
  
  mt <- tileMethylCounts(methylRawList.obj)
  expect_is(unite(mt, destrand = TRUE), "methylBase")
  
  # Create methylDB object
  mdb <- makeMethylDB(methylRawList.obj, dbdir = "methylDB")

  # Clean up
  on.exit(unlink("methylDB",recursive = TRUE), add = TRUE, after = TRUE)
  
  # Test base resolution with destranding
  expect_is(unite(mdb, destrand = TRUE), "methylBaseDB")
  
  # Create tiled version and set resolution to base to mimic strand "*"
  mdbt <- tileMethylCounts(mdb)
  # Test normal unite works for regions
  expect_is(unite(mdbt), "methylBaseDB")
  
  # Test normal unite works for regions
  expect_is(unite(mdbt, destrand = TRUE),"methylBaseDB")
  
  # Check if original tiled methylRawDB still exist 
  # see https://github.com/al2na/methylKit/issues/270
  expect_is(mdbt, "methylRawListDB")
  
  mdbt_base <- methylRawListDB(
    lapply(mdbt, function(x) {x@resolution <- "base"; x}),
    treatment = getTreatment(mdb)
  )
  
  # Test normal unite works
  expect_is(unite(mdbt_base), "methylBaseDB")
  
  # Test that no error is thrown when using destrand=TRUE with region resolution
  # see  https://github.com/al2na/methylKit/issues/271
  expect_is(unite(mdbt_base, destrand = TRUE),"methylBaseDB")
  
})

