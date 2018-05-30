context("SummarizedExperiment")
test_that("SummarizedExperiment can be created", {
  skip_if_devel()
  closeAll(sl)
  sumex <- as(sl[seq_len(978), seq_len(10)], "SummarizedExperiment")
  expect_equal(sumex$inst_id[1],  "CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17")
  expect_true(all.equal(SummarizedExperiment::assays(sumex)[[1]][5,10],  
                        12, 
                        tol = 0.00001))
  expect_equal(as.numeric(nrow(sumex)),  978)
  expect_equal(as.numeric(ncol(sumex)),  10)
  closeAll(sl)
})
test_that("loadL1K can load dataset without arguments", {
  skip_if_devel()
  closeAll(sl)
  sumex <- loadL1K(sl[,seq_len(20)])
  expect_equal(as.numeric(ncol(sumex)),  20)
  closeAll(sl)
})
test_that("SummarizedExperiment can be created by where clause", {
  skip_if_devel()
  where_clause = list("distil_id" = 'CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17')
  sumex <- loadL1K(sl, where_clause = where_clause)
  expect_equal(sumex$inst_id[1],  "CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17")
  expect_equal(as.numeric(nrow(sumex)),  12328)
  expect_equal(as.numeric(ncol(sumex)),  5)
  closeAll(sl)
})

test_that("SummarizedExperiment can be created with controls id'd automatically", {
  skip_if_devel()
  closeAll(sl)
  where_clause = list("distil_id" = 'CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17')
  sumex <- loadL1K(sl, where_clause = where_clause, controls = TRUE)
  expect_equal(sumex$inst_id[1],  "CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17")
  expect_equal(as.numeric(nrow(sumex)),  12328)
  expect_equal(as.numeric(ncol(sumex)),  64)
  closeAll(sl)
})

test_that("Failed query exits somewhat gracefully", {
  skip_if_devel()
  closeAll(sl)
  expect_warning( 
    tt <- loadL1K(sl, where_clause = list(pert_iname = "foobar")), 
    "no results")
  expect_null(tt)
  closeAll(sl)
})
