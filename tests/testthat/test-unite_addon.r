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


# does unite work as expected ?
samp1 <- data.frame(chr = "chr11", 
                    start = c(119000000,119000001),
                    end = c(119000000,119000001),
                    strand = c("+","-"),
                    coverage = c(17,14),
                    numCs = c(13,14),
                    numTs = c(4,0))

samp1Raw <- new("methylRaw",samp1,sample.id="samp1",assembly="none",
    context="none",resolution="base")

samp2 <- data.frame(chr = "chr11", 
                    start = c(119000000,119000001),
                    end = c(119000000,119000001),
                    strand = c("+","-"),
                    coverage = c(23,12),
                    numCs = c(20,11),
                    numTs = c(3,1))

samp2Raw <- new("methylRaw",samp2,sample.id="samp2",assembly="none",
    context="none",resolution="base")

res <- data.frame(chr = as.character("chr11"), 
                  start = as.integer(119000000),
                  end = as.integer(119000000),
                  strand = as.character("+"),
                  coverage1 = c(31),
                  numCs1 = c(27),
                  numTs1 = c(4),
                  coverage2 = c(35),
                  numCs2 = c(31),
                  numTs2 = c(4),
                  stringsAsFactors = FALSE
)


test_that("check if destranding works as expected", {
  expect_equal(getData(unite(
    methylRawList(samp1Raw, samp2Raw, treatment = c(0, 1)), destrand = TRUE
  )),
  res)
})


# ## benchmark for .CpG.dinuc.unify vs .CpG.dinuc.unifyOld
# df <- data.table::rbindlist(lapply(methylRawList.obj, getData))
# dim(df)
# # df <- rbind(df,df,df,df,df,df,df,df)
# # df$chr <- rep(x = 1:8,each = dd[1])
# for(i in 1:4 ) {
#   df2 <- copy(df)
#   df2[,`:=`(start = start - i,end = end - i)]
#   df <- rbind(df,df2)
# }
# data.table::setkey(df,chr,start,end)
# df <- df[!duplicated(df[,1:3]),]
# df <- as.data.frame(df)
# dim(df)
# 
# 
# bench::mark(
#   .CpG.dinuc.unify(df),
#   .CpG.dinuc.unifyOld(df)
# )

unlink(dbdir,recursive = TRUE)
