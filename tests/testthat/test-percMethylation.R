context("test percMethylation function")

data("methylKit")

dbdir = "methylDB"
methylBaseDB = makeMethylDB(methylBase.obj, dbdir = dbdir)

## no save txt
test_that("check if percMethylation is same for methylBase and methylBaseDB",
          {
            expect_equal(
              percMethylation(methylBase.obj, save.txt = FALSE, rowids = TRUE),
              percMethylation(methylBaseDB, save.txt = FALSE, rowids = TRUE)
            )
            expect_equal(
              percMethylation(methylBase.obj, save.txt = FALSE, rowids = FALSE),
              percMethylation(methylBaseDB, save.txt = FALSE, rowids = FALSE)
            )
          })

## save txt
test_that("check if percMethylation is writting to file for methylBaseDB",
          {
            expect_equal(
              as.matrix(read.table(
                percMethylation(methylBaseDB, save.txt = TRUE, rowids = TRUE),
                header = TRUE,
                row.names = 1
              )),
              percMethylation(methylBase.obj, save.txt = TRUE, rowids = TRUE)
            )
            
            expect_equal(
              as.matrix(read.table(
                percMethylation(methylBaseDB, save.txt = TRUE, rowids = FALSE),
                header = TRUE
              )),
              percMethylation(methylBase.obj, save.txt = TRUE, rowids = FALSE)
            )
          })


unlink(dbdir, recursive = TRUE)