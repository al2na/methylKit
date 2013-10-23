context("unite checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list,
            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
            pipeline="amp",treatment=c(1,1,0,0))

# unite function
ta=unite(myobj,destrand=T)
methidh=unite(myobj)
methidh2=unite(myobj,min.per.group=1L)


test_that("test if output of unite is  methylBase object", {
    expect_that(ta, is_a('methylBase'))
    expect_that(methidh, is_a('methylBase'))
    expect_that(methidh2, is_a('methylBase'))

})
