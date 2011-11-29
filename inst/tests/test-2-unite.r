context("unite checks")

file.list=list( system.file("data", "test1.myCpG.txt", package = "methylKit"),
                system.file("data", "test2.myCpG.txt", package = "methylKit"),
                system.file("data", "control1.myCpG.txt", package = "methylKit"),
                system.file("data", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list,
            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

# unite function
t=unite(myobj,destrand=T)
methidh=unite(myobj)


test_that("test if output of unite is  methylBase object", {
    expect_that(t, is_a('methylBase'))
    expect_that(methidh, is_a('methylBase'))
})
