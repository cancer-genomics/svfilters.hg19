context("Germline filters")
test_that("germlineFilters", {
  gfs <- listGenomeFilters("hg19")
  expect_is(gfs, "list")
  expect_error(listGenomeFilters("hg18"))
})
