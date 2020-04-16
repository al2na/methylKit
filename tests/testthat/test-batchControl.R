context("check if batch control functions work")

data(methylKit)

dbdir <- "methylDB"

methylBaseDB = makeMethylDB(methylBase.obj, dbdir = dbdir)
sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),age=c(19,34,23,40))

mat=percMethylation(methylBase.obj)
mat[mat==100]=80

matFail = percMethylation(methylBase.obj)
matFail[matFail==max(matFail)]=200

matFail2 = percMethylation(methylBase.obj)
matFail2 = cbind(matFail2,1)


test_that("check if assocComp works", {
  expect_length(assocComp(mBase = methylBase.obj, sampleAnnotation), 3)
  expect_error(
    assocComp(mBase = methylBase.obj, data.frame()),
    "Sample Annotation has to have at least one column."
  )
  expect_equal(
    assocComp(mBase = methylBase.obj, sampleAnnotation),
    assocComp(mBase = methylBaseDB, sampleAnnotation)
  )
})


test_that("check if removeComp works", {
  expect_is(removeComp(methylBase.obj, comp = 1), "methylBase")
  expect_is(removeComp(methylBaseDB, comp = 1), "methylBaseDB")
  expect_is(removeComp(methylBase.obj, comp = c(3, 4)), "methylBase")
  expect_is(removeComp(methylBaseDB, comp = c(3, 4)), "methylBaseDB")
  expect_error(removeComp(methylBase.obj), "no component to remove")
  expect_error(removeComp(methylBaseDB), "no component to remove")
  expect_error(
    removeComp(methylBase.obj, comp = 10),
    "elements can only take values between 1 and number of samples"
  )
  expect_error(
    removeComp(methylBase.obj, comp = 10),
    "elements can only take values between 1 and number of samples"
  )
  
})

test_that("check if reconstruct works", 
  {
    expect_is(reconstruct(mat, methylBase.obj), "methylBase")
    expect_is(reconstruct(mat, methylBaseDB), "methylBaseDB")
    expect_is(reconstruct(mat, methylBaseDB, save.db = TRUE, dbdir = dbdir),
              "methylBaseDB")
    expect_is(reconstruct(mat, methylBaseDB, save.db = FALSE), "methylBase")
    expect_error(
      reconstruct(matFail, methylBase.obj),
      "make sure 'methMat' is percent methylation matrix"
    )
    expect_error(
      reconstruct(matFail, methylBaseDB),
      "make sure 'methMat' is percent methylation matrix"
    )
    expect_error(
      reconstruct(matFail2, methylBase.obj),
      "methMat dimensions do not match number of samples"
    )
    expect_error(
      reconstruct(matFail2, methylBaseDB),
      "methMat dimensions do not match number of samples"
    )
    
})

unlink(dbdir,recursive = TRUE)