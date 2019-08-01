context("test reading of text and tabix files")

# I have two files, one is already processed, the other not
compr.file.path <- system.file("extdata", "test1.txt.bgz", package = "methylKit")
uncompr.file.path <- system.file("extdata","test1.myCpG.txt",package = "methylKit")

dbdir="methylDB"

# the compressed can be directly loaded by using the path to the database file
test_that("reading of tabix without dbtype leads to error", {
  expect_error(methRead(location = compr.file.path,sample.id = "test1",assembly = "hg18",dbdir = dbdir))
})

# but you need to take care of the dbtype parameter ( is NA by default)
myobj.db <- methRead(location = compr.file.path,sample.id = "test1",assembly = "hg18",dbtype = "tabix",dbdir = dbdir)
test_that("if read of database return a methylRawDB", {
  expect_is(myobj.db, 'methylRawDB')
})


# the uncompressed needs to be processed before
myobj2.db <- methRead(location = uncompr.file.path,sample.id = "test1",assembly = "hg18",dbtype = "tabix",dbdir = dbdir)
test_that("if read of database return a methylRawDB", {
  expect_is(myobj2.db, 'methylRawDB')
})


# merging uncompressed and compressed files will lead to error
test_that("reading of merged uncompressed and compressed files leads to error", {
  expect_error(methRead(location = list(compr.file.path,uncompr.file.path),sample.id = list("test1","test2"),treatment = c(1,1),assembly = "hg18",dbtype = "tabix",dbdir = dbdir))
})


# you have to process uncompressed files beforehand and can than load them together
myobj2.list.db <- methRead(location = list(compr.file.path,myobj2.db@dbpath),sample.id = list("test1","test2"),treatment = c(1,1),assembly = "hg18",dbtype = "tabix",dbdir = dbdir)
test_that("reading of compressed files returns methylRawListDB", {
  expect_is(myobj2.list.db, 'methylRawListDB')
})


unlink("tests/testthat/methylDB",recursive = TRUE)