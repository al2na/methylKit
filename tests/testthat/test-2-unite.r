context("unite checks")

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
ta=unite(myobj,destrand=T)
methidh=unite(myobj)
methidh2=unite(myobj,min.per.group=1L)

suppressMessages(taDB <- unite(mydblist,destrand=T))
suppressMessages(methidhDB <- unite(mydblist))
suppressMessages(methidh2DB <- unite(mydblist,min.per.group=1L))


suppressMessages(ta2methylDB <- unite(myobj,destrand=T,save.db = T,dbdir="methylDB"))
taDB2methyl <- unite(mydblist,destrand=T,save.db = F)


test_that("test if output of unite is  methylBase object", {
    expect_is(ta, 'methylBase')
    expect_is(methidh, 'methylBase')
    expect_is(methidh2, 'methylBase')
    expect_is(taDB2methyl, 'methylBase')

})

test_that("test if output of unite is  methylBaseDB object", {
  expect_is(taDB, 'methylBaseDB')
  expect_is(methidhDB, 'methylBaseDB')
  expect_is(methidh2DB, 'methylBaseDB')
  expect_is(ta2methylDB, 'methylBaseDB')
})

