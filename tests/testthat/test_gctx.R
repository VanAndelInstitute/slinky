context("GCTX file I/O")
test_that("Column names can be retrieved", {
  skip_if_devel()
  tt <- colnames(sl[,seq_len(3)])
  expect_equal(length(tt),  3)
  expect_equal(tt[1], "CPC020_MCF7_24H_X1_F1B4_DUO52HI53LO:P17")
})
test_that("Row names can be retrieved", {
  skip_if_devel()
  tt <- rownames(sl[seq_len(3),])
  expect_equal(length(tt),  3)
  expect_equal(tt[1], "5720")
})
test_that("Matrix data can be read", {
  skip_if_devel()
  tt <- readGCTX(sl[seq_len(3), seq_len(6) + 4])
  expect_equal(nrow(tt), 3)
  expect_equal(ncol(tt), 6)
  expect_equal(colnames(tt)[6], "KDA001_MCF7_144H_X2_B1_DUO44HI45LO:P15")
})
test_that("Correct columns are selected", {
  skip_if_devel()
  # load data directly from datafile
  ids.treat <- colnames(sl[, seq_len(5)])
  ids.ctrl <- colnames(sl)[seq_len(59) + 17]
  data <- rhdf5::h5read(gctx, name = "0/DATA/0/matrix")
  coln <- rhdf5::h5read(gctx, name = "0/META/COL/id" )
  data.ctrl <- data[, match(ids.ctrl, coln)]
  data.treat <- data[, match(ids.treat, coln)]

  # compare to columns extracted by colnames
  ix <- match(c(ids.treat, ids.ctrl), colnames(sl))
  expect_equal(as.character(c(ids.treat, ids.ctrl)), 
          as.character(colnames(sl)[ix]))
  closeAll(sl)
})
