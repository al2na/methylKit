context("clusterSamples and PCASamples checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(read( file.list,
                                  sample.id=list("t1","t2","c1","c2"),assembly="hg18",
                                  pipeline="amp",treatment=c(1,1,0,0),dbtype = "tabix",dbdir="methylDB"))

# unite function
methidh=unite(myobj)
suppressMessages(methidhDB <- unite(mydblist))

test_that("check if clusterSamples output is a hclust tree object", {
    expect_is(clusterSamples(methidh, dist="correlation", method="ward", plot=FALSE), 
        'hclust')
    expect_is(clusterSamples(methidhDB, dist="correlation", method="ward", plot=FALSE), 
            'hclust')
})

test_that("check if PCASamples output is a summary of princomp", {
    expect_is(PCASamples(methidh,obj.return=TRUE),
        'prcomp')
    expect_is(PCASamples(methidhDB,obj.return=TRUE),
            'prcomp')
})
