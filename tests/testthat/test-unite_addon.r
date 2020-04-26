context("check if some unite works")

data(methylKit)
dbdir <- "methylDB"


methylRawListDB = makeMethylDB(methylRawList.obj, dbdir = dbdir)

test_that("check if unite works", {
  expect_equal(
    getData(unite(methylRawList.obj)),
    getData(unite(methylRawListDB))
  )
  expect_equal(
    getData(unite(methylRawList.obj,destrand = TRUE)),
    getData(unite(methylRawListDB,destrand = TRUE))
  )
  expect_equal(
    getData(unite(methylRawList.obj,destrand = FALSE,min.per.group = 1L)),
    getData(unite(methylRawListDB,destrand = FALSE,min.per.group = 1L))
  )
  expect_equal(
    getData(unite(methylRawList.obj,destrand = TRUE,min.per.group = 1L)),
    getData(unite(methylRawListDB,destrand = TRUE,min.per.group = 1L))
  )
})

## benchmark for .CpG.dinuc.unify vs .CpG.dinuc.unifyOld

# df <- rbindlist(lapply(methylRawList.obj,getData))
# setkey(df,chr,start,end)
# df <- df[!duplicated(df[,1:3]),]
# df <- as.data.frame(df)
# 
bench::mark(
  .CpG.dinuc.unify(df),
  .CpG.dinuc.unifyOld(df)
)

microbenchmark::microbenchmark(
  .CpG.dinuc.unify(df),
  .CpG.dinuc.unifyOld(df)
)


unlink(dbdir,recursive = TRUE)
