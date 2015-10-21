
library(methylKit)
file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

methylRawList.obj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",treatment=c(1,1,0,0))


methylBase.obj=unite(methylRawList.obj)

methylDiff.obj=calculateDiffMeth(methylBase.obj)

db.list=list( system.file("extdata", "test1.txt.bgz", package = "methylKit"),
                system.file("extdata", "test2.txt.bgz", package = "methylKit"),
                system.file("extdata", "ctrl1.txt.bgz", package = "methylKit"),
                system.file("extdata", "ctrl2.txt.bgz", package = "methylKit") )

methylRawListDB.obj=read( db.list,
                          sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",treatment=c(1,1,0,0),
                          dbtype = "tabix")

methylBaseDB.obj=unite(methylRawListDB.obj)

save(methylRawList.obj,methylBase.obj,methylDiff.obj,methylRawListDB.obj,methylBaseDB.obj,file="methylKit.RData")
