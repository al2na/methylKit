context("tests related to tileMethylCounts")



data(methylKit)
output_path_issue <- file.path(tempdir(),"special_folder") |> normalizePath()
getwd()

# simulated https://github.com/al2na/methylKit/issues/266
tiles_raw.issue <- tileMethylCounts(
  object = methylRawList.obj[[1]],
  win.size = 500,
  step.size = 500,
  save.db = TRUE,
  suffix = "tiled",
  dbdir = file.path(output_path_issue)
) 

tiles_rawlist.issue <- tileMethylCounts(
  object = methylRawList.obj,
  win.size = 500,
  step.size = 500,
  save.db = TRUE,
  suffix = "tiled",
  dbdir = file.path(output_path_issue)
) 

tiles_base.issue <- tileMethylCounts(
  object = methylBase.obj,
  win.size = 500,
  step.size = 500,
  save.db = TRUE,
  suffix = "tiled",
  dbdir = file.path(output_path_issue)
)

tiles_rawdb.issue <- tileMethylCounts(
  object = makeMethylDB(dbdir = tempdir(), methylRawList.obj[[1]]),
  win.size = 500,
  step.size = 500,
  save.db = TRUE,
  suffix = "tiled",
  dbdir = file.path(output_path_issue)
) 


tiles_rawlistdb.issue <- tileMethylCounts(
  object = makeMethylDB(dbdir = tempdir(), methylRawList.obj),
  win.size = 500,
  step.size = 500,
  save.db = TRUE,
  suffix = "tiled",
  dbdir = file.path(output_path_issue)
) 

tiles_basedb.issue <- tileMethylCounts(
  object = makeMethylDB(dbdir = tempdir(), methylBase.obj),
  win.size = 500,
  step.size = 500,
  save.db = TRUE,
  suffix = "tiled",
  dbdir = file.path(output_path_issue)
)

test_that("test if output of tiling methylRawList is saved in correct path ", {
  expect_identical(
    tiles_rawlist.issue |> getDBPath() |> normalizePath() |> dirname() |> unique(),
    output_path_issue |> normalizePath()
  )
})
test_that("test if output of tiling methylRaw is saved in correct path ", {
  expect_identical(
    tiles_raw.issue |> getDBPath() |> normalizePath() |>  dirname()  |> unique(),
    output_path_issue |> normalizePath()
  )
})
test_that("test if output of tiling methylBase is saved in correct path ", {
  expect_identical(
    tiles_base.issue |> getDBPath() |> normalizePath() |> dirname()|> unique(),
    output_path_issue |> normalizePath()
  )
})

test_that("test if output of tiling methylRawList is saved in correct path ", {
  expect_identical(
    tiles_rawlistdb.issue |> getDBPath() |> normalizePath() |> dirname() |> unique(),
    output_path_issue |> normalizePath()
  )
})
test_that("test if output of tiling methylRaw is saved in correct path ", {
  expect_identical(
    tiles_rawdb.issue |> getDBPath() |> normalizePath() |>  dirname()  |> unique(),
    output_path_issue |> normalizePath()
  )
})
test_that("test if output of tiling methylBase is saved in correct path ", {
  expect_identical(
    tiles_basedb.issue |> getDBPath() |> normalizePath() |> dirname()|> unique(),
    output_path_issue |> normalizePath()
  )
})
