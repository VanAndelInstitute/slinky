# usage:  testthat::auto_test("R/", "tests/testthat/")
# or ... test_dir("tests/testthat")

# fetch the clue demo key for testing purposes.  Please use your own key for any production work.
suppressMessages(library(Biobase))
user_key <- httr::content(httr::GET("https://api.clue.io/temp_api_key"), as="parsed")$user_key
sl <- Slinky$new(user_key)

context("Class instantiation")
test_that("Object can be created", {  
  # don't access private slots in real life
  expect_is(sl,  "Slinky")
})

context("Count function")
test_that("Record count can be retrieved", {  
  tt <- sl$clue.count("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(tt, 1)
})
test_that("Records can be fetched", {  
  tt <- sl$clue("genes", fields=c("entrez_id", "gene_symbol"), limit=1)
  expect_equal(nrow(tt), 1)
})

context("Genes API")
test_that("A gene can be retrieved", {  
  tt <- sl$clue("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(nrow(tt),  1)
})
test_that("Limit working correctly on genes API", {  
  tt <- sl$clue("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(nrow(tt), 1)
  tt <- sl$clue("genes", fields=c("entrez_id", "gene_symbol", "gene_name"), limit=20)
  expect_equal(nrow(tt), 20)
})
test_that("Many genes can be retrieved using a cluster object", {  
  skip("Skipping clustered test because check environment does not like it.")
  skip_on_cran()
  cl <- parallel::makeCluster(4)
  tt <- sl$clue("genes", fields=c("entrez_id", "gene_symbol", "gene_name"), limit=2000, cl=cl)
  expect_equal(nrow(tt),  2000)
  expect_equal(tt[1,2], "FAM200B")
  parallel::stopCluster(cl)
})

context("Cells API")
test_that("Cell lines can be retrieved by iname", {  
  tt <- sl$clue("cells", fields=c("cell_iname", "cell_lineage"), 
                 where_clause=list("cell_id"=list(like="^A3")))
  expect_equal(nrow(tt),  3)
})
  
context("Signatures API")
test_that("A signature can be retrieved", {  
  tt <- sl$clue("sigs", where_clause=list(pert_desc="sirolimus"), fields = c("is_gold", "distil_id"), unpack_sigs = FALSE)
  expect_equal(nrow(tt),  202)
})
test_that("Instances can be unpacked from signature endpoint", {  
  tt <- sl$clue("sigs", where_clause=list(pert_desc = "sirolimus"), fields = c("is_gold", "distil_id"), unpack_sigs = TRUE)
  expect_equal(nrow(tt),  664)
})
test_that("A list of instance ids can be retrieved", {
  tt <- sl$clue.instances(where_clause=list("pert_desc" = "sirolimus", "is_gold" = TRUE, cell_id = 'MCF7'))
  expect_equal(length(tt), 62)
})
test_that("Plate ids can be extracted from dataframe", {
  tt <- sl$clue("sigs", where_clause=list(pert_desc = "sirolimus"), fields = c("is_gold", "distil_id"), unpack_sigs = TRUE)
  tt <- sl$distil.to.plate(tt)
  expect_equal(length(tt), 664)
  expect_equal(tt[1], "ASG001_MCF7_24H_X1")
})

context("Eset")
test_that("Eset can be created", {  
  sl$close()
  sl <- Slinky$new(user_key, system.file("extdata", "demo.gctx", package="slinky"))
  eset <- sl$toEset(index = list(1:978, 1:10), info_file = system.file("extdata", "demo_inst_info.txt", package="slinky"))
  expect_equal(pData(eset)[1,1],  "CPC004_A375_6H_X1_B3_DUO52HI53LO:H04")
  expect_true(all.equal(exprs(eset)[5,10],  7.7316, tol=0.00001))
  expect_equal(as.numeric(nrow(eset)),  978)
  expect_equal(as.numeric(ncol(eset)),  10)
  sl$close()
})

context("GCTX file I/O")
gctx <- system.file("extdata", "demo.gctx", package="slinky")
sl <- Slinky$new(user_key, gctx)
test_that("Column names can be retrieved", {  
  tt <- sl$colnames(index=list(1:3))
  expect_equal(length(tt),  3)
  expect_equal(tt[1], "CPC004_A375_6H_X1_B3_DUO52HI53LO:H04")
})
test_that("Row names can be retrieved", {  
  tt <- sl$rownames(index=list(1:3))
  expect_equal(length(tt),  3)
  expect_equal(tt[1], 5720)
})
test_that("Matrix data can be read", {  
  tt <- sl$readGCTX(index=list(1:3, 5:10))
  expect_equal(nrow(tt), 3)
  expect_equal(ncol(tt), 6)
  expect_equal(colnames(tt)[6], "CPC004_A375_6H_X1_B3_DUO52HI53LO:H20")
})
sl$close()
test_that("Matrix data can be read with file specified by constructor", {  
  sl = Slinky$new(key = user_key, gctx = system.file("extdata", "demo.gctx", package="slinky"))
  tt <- sl$readGCTX(gctx, index=list(1:3, 5:10))
  expect_equal(nrow(tt), 3)
  expect_equal(ncol(tt), 6)
  expect_equal(colnames(tt)[6], "CPC004_A375_6H_X1_B3_DUO52HI53LO:H20")
})

context("Characteristic direction")
sl <- Slinky$new(user_key, system.file("extdata", "demo.gctx", package="slinky"))
eset <- sl$toEset(index = list(1:978, 1:10), info_file = system.file("extdata", "demo_inst_info.txt", package="slinky"))
test_that("Characteristic direction can be calculated based on two expression sets", {  
  tt <- sl$chDir(eset[, 1:4], eset[, 5:10])
  expect_equivalent(tt[1], 0.017179161132863912664348)
})
test_that("Characteristic direction can be calculated based on two matrices", {  
  tt <- sl$chDir(exprs(eset[, 1:4]), exprs(eset[, 5:10]))
  expect_equivalent(tt[1], 0.017179161132863912664348)
})


context("Enrichment score calc")
zsvc <- system.file("extdata", "test_inst_zsvc_x1000.rds", package="slinky")
info <- system.file("extdata", "test_inst_info.rds", package="slinky")
test_that("KS statistic can be calculated", {  
  data <- readRDS(zsvc) / 1000
  drugs <- readRDS(info)
  up <- rownames(data)[order(data[,1], decreasing = TRUE)][1:25]
  down <- rownames(data)[order(data[,1])][1:25]
  s <- sl$score(data, "ks", up=up, down=down)
  expect_equivalent(s[1], 1)
  # expect_equivalent(s[2],  0.47485834207765)
})
test_that("xsum statistic can be calculated", {  
  data <- readRDS(zsvc) / 1000
  drugs <- readRDS(info)
  up <- rownames(data)[order(data[,1], decreasing = TRUE)][1:25]
  down <- rownames(data)[order(data[,1])][1:25]
  s <- sl$score(data, "xsum", up=up, down=down)
  expect_equivalent(s[1], 185.634)
})
