context("test file list and modRead and getMethylationStats check")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=modRead( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

test_that("check if there are 4 test files in the file.list",{
    expect_that(length(file.list),
        equals(4))
})

test_that("if read return a methylRawlist", {
    expect_that(myobj, is_a('methylRawList'))
})

test_that("getMethylationStats on myobj[[2]] works", {
    expect_that (getMethylationStats(myobj[[2]],plot=F,both.strands=F),
            prints_text('methylation statistics per base'))
})

