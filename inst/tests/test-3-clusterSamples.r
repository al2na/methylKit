context("clusterSamples and PCASamples checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

# unite function
methidh=unite(myobj)

test_that("check if clusterSamples output is a hclust tree object", {
    expect_that(clusterSamples(methidh, dist="correlation", method="ward", plot=FALSE), 
        is_a('hclust'))
})

test_that("check if PCASamples output is a summary of princomp", {
    expect_that(PCASamples(methidh,obj.return=TRUE),
        is_a('prcomp'))
})
