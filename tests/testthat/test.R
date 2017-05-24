# usage:  testthat::auto_test("R/", "tests/testthat/")
# or ... test_dir("tests/testthat")

library(doParallel)
Sys.setenv(CLUE_API_KEY="db1db72b80bcb8ed812a2a7af0198bbb")

sl <- Slinky$new()

context("Class instantiation")
test_that("Object can be created", {  
  # don't access private slots in real life
  expect_is(sl,  "Slinky")
})

context("Count function")
test_that("Record count can be retrieved", {  
  tt <- sl$count("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(tt, 1)
})
test_that("Records can be fetched", {  
  tt <- sl$fetch("genes", fields=c("entrez_id", "gene_symbol"), limit=1)
  expect_equal(nrow(tt), 1)
})

context("Genes API")
test_that("A gene can be retrieved", {  
  tt <- sl$fetch("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(nrow(tt),  1)
})
test_that("Limit working correctly on genes API", {  
  tt <- sl$fetch("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(nrow(tt), 1)
  tt <- sl$fetch("genes", fields=c("entrez_id", "gene_symbol", "gene_name"), limit=20)
  expect_equal(nrow(tt), 20)
})
test_that("Many genes can be retrieved using a cluster object", {  
  cl <- makeCluster(4)
  tt <- sl$fetch("genes", fields=c("entrez_id", "gene_symbol", "gene_name"), limit=2000, cl=cl)
  expect_equal(nrow(tt),  2000)
  expect_equal(tt[1,2], "FAM200B")
  stopCluster(cl)
})

context("Cells API")
test_that("Cell lines can be retrieved by iname", {  
  tt <- sl$fetch("cells", fields=c("cell_iname", "cell_lineage"), 
                 where_clause=list("cell_id"=list(like="^A3")))
  expect_equal(nrow(tt),  3)
})
  
context("Signatures API")
test_that("A signature can be retrieved", {  
  tt <- sl$fetch("sigs", where_clause=list(pert_desc="sirolimus"), fields=c("is_gold", "distil_id"))
  expect_equal(nrow(tt),  202)
})
