context("test methReadList and methReadListDB constructors")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

list_myobj = mapply(methRead,
                    location=file.list,
                    sample.id=c("test1","test2","ctrl1","ctrl2"),
                    MoreArgs = list(assembly="hg18", pipeline="amp"),
                    SIMPLIFY = FALSE
                    )

myobjlist_joined <- methylRawList(list_myobj, treatment=c(1,1,0,0))

myobjlist=methRead( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

list_mydb = suppressMessages(mapply(methRead,
                   location=file.list,
                   sample.id=c("test1","test2","ctrl1","ctrl2"),
                   MoreArgs = list(assembly="hg18",pipeline="amp",
                   dbtype = "tabix",dbdir="methylDB"),
                   SIMPLIFY = FALSE))

mydblist_joined <- methylRawListDB(list_mydb,treatment=c(1,1,0,0))

mydblist = suppressMessages(methRead( file.list,
                                      sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0),dbtype = "tabix",dbdir="methylDB"))



test_that("check if the joined list is equal to direct processing for methylRawList",{
  expect_identical(myobjlist_joined,myobjlist)
})

test_that("check if the joined list is equal to direct processing for methylRawListDB",{
  expect_identical(mydblist_joined,mydblist)
})


unlink("tests/testthat/methylDB",recursive = TRUE)