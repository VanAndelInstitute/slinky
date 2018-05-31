context("Count function")
test_that("Record count can be retrieved", {
  skip_if_devel()
  tt <- clueCount(sl, "genes", where_clause = list(gene_symbol = "TP53"))
  expect_equal(tt, 1)
})
test_that("Records can be fetched", {
  skip_if_devel()
  tt <- clue(sl, "genes", fields = c("entrez_id", "gene_symbol"), limit = 1)
  expect_equal(nrow(tt), 1)
})

context("Genes API")
test_that("A gene can be retrieved", {
  skip_if_devel()
  tt <- clue(sl, "genes", where_clause = list(gene_symbol = "TP53"))
  expect_equal(nrow(tt),  1)
})
test_that("Limit working correctly on genes API", {
  skip_if_devel()
  tt <- clue(sl, "genes", where_clause = list(gene_symbol = "TP53"))
  expect_equal(nrow(tt), 1)
  tt <- clue(sl, "genes", fields = c("entrez_id", "gene_symbol",
                                  "gene_name"), limit = 20)
  expect_equal(nrow(tt), 20)
})

context("Cells API")
test_that("Cell lines can be retrieved by iname", {
  skip_if_devel()
  tt <- clue(sl, "cells", fields = c("cell_iname", "cell_lineage"),
                 where_clause = list("cell_id" = list(like = "^A3")))
  expect_equal(nrow(tt),  3)
})

context("Signatures API")
test_that("A signature can be retrieved", {
  skip_if_devel()
  tt <- clue(sl, "sigs", where_clause = list(pert_desc = "sirolimus"),
                fields = c("is_gold", "distil_id"),
                unpack_sigs = FALSE)
  expect_equal(nrow(tt),  202)
})
test_that("Instances can be unpacked from signature endpoint", {
  skip_if_devel()
  tt <- clue(sl, "sigs", where_clause = list(pert_desc = "sirolimus"),
                fields = c("is_gold", "distil_id"),
                unpack_sigs = TRUE)
  expect_equal(nrow(tt),  664)
})
test_that("A list of instance ids can be retrieved", {
  skip_if_devel()
  tt <- clueInstances(sl, where_clause = list("pert_desc" = "sirolimus",
                                            "is_gold" = TRUE,
                                            cell_id = 'MCF7'))
  expect_equal(length(tt), 62)
})

context("Profiles API")
test_that("Correct controls are selected", {
  skip_if_devel()
  # load data directly from datafile
  ids.ctrl <- colnames(sl)[18:76]
  data <- rhdf5::h5read(gctx, name = "0/DATA/0/matrix")
  coln <- rhdf5::h5read(gctx, name = "0/META/COL/id" )
  data.ctrl <- data[, match(ids.ctrl, coln)]
  ids.ctrl2 <- controls(sl, 
                        metadata(sl[,seq_len(5)])$distil_id)$distil_id
  expect_equal(as.character(ids.ctrl2), as.character(ids.ctrl))
})
test_that("Instance ids can be retrieved", {
  skip_if_devel()
  tt <- clueInstances(sl, where_clause = list(pert_desc = "sirolimus"))
  expect_equal(length(tt), 503)
})

