context("bins")
test_that("bins", {
  data(bins1kb_hg19)
  expect_identical(length(bins1kb_hg19),
                   2691009L)

  data(bins1kb_hg18)
  expect_identical(length(bins1kb_hg18),
                   2650192L)
})
