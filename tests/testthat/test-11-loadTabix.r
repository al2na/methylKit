context("test loading of tabix files")


file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list,
            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
            pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(read( file.list,
                                  sample.id=list("t1","t2","c1","c2"),assembly="hg18",
                                  pipeline="amp",treatment=c(1,1,0,0),dbtype = "tabix",dbdir="methylDB"))

# unite function
methidh=unite(myobj)

methdiff <- calculateDiffMeth(methidh)

write2tabix(myobj[[1]],file.name = "methylDB/my_raw.txt",rm.txt = FALSE)
write2tabix(mydblist[[1]],file.name = "methylDB/my_raw2.txt",rm.txt = FALSE)
write2tabix(methidh,file.name = "methylDB/my_base.txt",rm.txt = FALSE)
write2tabix(methdiff,file.name = "methylDB/my_diff.txt",rm.txt = FALSE)

no_header <- system.file("extdata", "ctrl1.txt.bgz", package = "methylKit")

# the compressed can be directly loaded by using the path to the database file
test_that("reading of tabix without heading leads to error", {
  expect_error(loadMethylRawDB(dbpath = no_header))
})

# the compressed can be directly loaded by using the path to the database file
raw <- loadMethylRawDB(dbpath =  "methylDB/my_raw.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(raw,'methylRawDB')
})

# the compressed can be directly loaded by using the path to the database file
raw2 <- loadMethylRawDB(dbpath =  "methylDB/my_raw2.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(raw2,'methylRawDB')
})

# the compressed can be directly loaded by using the path to the database file
base <- loadMethylBaseDB(dbpath =  "methylDB/my_base.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(base,'methylBaseDB')
})

# the compressed can be directly loaded by using the path to the database file
diff <- loadMethylDiffDB(dbpath =  "methylDB/my_diff.txt.bgz")
test_that("reading of tabix without dbtype leads to error", {
  expect_is(diff,'methylDiffDB')
})

