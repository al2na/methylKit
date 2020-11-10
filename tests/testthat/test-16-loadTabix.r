context("test loading of tabix files")


file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=methRead( file.list,
            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
            pipeline="amp",treatment=c(1,1,0,0))

suppressMessages({mydblist = methRead( file.list,
                     sample.id=list("t1","t2","c1","c2"),assembly="hg18",
                     pipeline="amp",treatment=c(1,1,0,0),
                     dbtype = "tabix",dbdir="methylDB")})

# process using usual functions
methidh=unite(myobj)
methidh_region <- tileMethylCounts(methidh)
methdiff <- calculateDiffMeth(methidh,save.db = TRUE)
methdiff_region <- calculateDiffMeth(methidh_region,save.db = TRUE)

methidh_db=unite(mydblist)
methidh_db_region <- tileMethylCounts(methidh_db)
methdiff_db <- calculateDiffMeth(methidh_db,save.db = TRUE)
methdiff_db_region <- calculateDiffMeth(methidh_db_region,save.db = TRUE)


no_header <- system.file("extdata", "ctrl1.txt.bgz", package = "methylKit")
# the compressed can be directly loaded by using the path to the database file
test_that("reading of tabix without heading leads to error", {
  expect_error(readMethylRawDB(dbpath = no_header))
})

# the compressed can be directly loaded by using the path to the database file
obj2tabix(myobj[[1]],filename = "methylDB/my_raw.txt",rm.txt = FALSE)
raw <- readMethylRawDB(dbpath =  "methylDB/my_raw.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(raw,'methylRawDB')
})

# the compressed can be directly loaded by using the path to the database file
obj2tabix(mydblist[[1]],filename = "methylDB/my_raw2.txt",rm.txt = FALSE)
raw2 <- readMethylRawDB(dbpath =  "methylDB/my_raw2.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(raw2,'methylRawDB')
})

# the compressed can be directly loaded by using the path to the database file
obj2tabix(methidh,filename = "methylDB/my_base.txt",rm.txt = FALSE)
base <- readMethylBaseDB(dbpath =  "methylDB/my_base.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(base,'methylBaseDB')
})

# the compressed can be directly loaded by using the path to the database file
obj2tabix(methidh_region,filename = "methylDB/my_base2.txt",rm.txt = FALSE)
base <- readMethylBaseDB(dbpath =  "methylDB/my_base2.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(base,'methylBaseDB')
})

# the compressed can be directly loaded by using the path to the database file
obj2tabix(methdiff,filename = "methylDB/my_diff.txt",rm.txt = FALSE)
diff <- readMethylDiffDB(dbpath =  "methylDB/my_diff.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(diff,'methylDiffDB')
})

test_that("reading of tabix can be done with one wrapper readMethylDB", {
  expect_error(readMethylDB(dbpath = no_header))
  expect_is(readMethylDB(dbpath =  "methylDB/my_raw.txt.bgz"),'methylRawDB')
  expect_is(readMethylDB(dbpath =  mydblist[[1]]@dbpath),'methylRawDB')
  expect_is(readMethylDB(dbpath =  "methylDB/my_base.txt.bgz"),'methylBaseDB')
  expect_is(readMethylDB(dbpath =  methidh_db@dbpath),'methylBaseDB')
  expect_is(readMethylDB(dbpath =  "methylDB/my_diff.txt.bgz"),'methylDiffDB')
  expect_is(readMethylDB(dbpath =  methdiff_db@dbpath),'methylDiffDB')
})


unlink("tests/testthat/methylDB",recursive = TRUE)
