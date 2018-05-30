context("Class instantiation")
test_that("Object can be created", {
  expect_is(sl,  "Slinky")
})
test_that("Dim functions working", {
  skip_if_devel()
  expect_is(sl,  "Slinky")
  expect_equal(ncol(sl), 131)
  expect_equal(nrow(sl), 12328)
  expect_equal(length(colnames(sl)), 131)
  expect_equal(length(rownames(sl)), 12328)
  expect_equal(rownames(sl)[1], "5720")
  expect_equal(colnames(sl)[1], "CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17")
})
