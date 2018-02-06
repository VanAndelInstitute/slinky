# don't rerun all tests when we are working on somethingin particular. 
# devel should be set to FALSE for deployment
devel <- FALSE

skip_if_devel <- function() {
  if (devel) {
    skip("Skipping to speed up development")
  }
}

# fetch the clue demo key so we can execuet our test suite.  
# Do NOT use the demo key for any production work.
suppressMessages(library(Biobase))
user_key <- httr::content(httr::GET("https://api.clue.io/temp_api_key"), 
                          as="parsed")$user_key
sl <- Slinky$new(user_key)

context("Class instantiation")
test_that("Object can be created", {  
  # don't access private slots in real life
  expect_is(sl,  "Slinky")
})

context("Count function")
test_that("Record count can be retrieved", {  
  skip_if_devel()
  tt <- sl$clue.count("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(tt, 1)
})
test_that("Records can be fetched", {  
  skip_if_devel()
  tt <- sl$clue("genes", fields=c("entrez_id", "gene_symbol"), limit=1)
  expect_equal(nrow(tt), 1)
})

context("Genes API")
test_that("A gene can be retrieved", {  
  skip_if_devel()
  tt <- sl$clue("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(nrow(tt),  1)
})
test_that("Limit working correctly on genes API", {  
  skip_if_devel()
  tt <- sl$clue("genes", where_clause=list(gene_symbol="TP53"))
  expect_equal(nrow(tt), 1)
  tt <- sl$clue("genes", fields=c("entrez_id", "gene_symbol", 
                                  "gene_name"), limit=20)
  expect_equal(nrow(tt), 20)
})

context("Cells API")
test_that("Cell lines can be retrieved by iname", {  
  skip_if_devel()
  tt <- sl$clue("cells", fields=c("cell_iname", "cell_lineage"), 
                 where_clause=list("cell_id"=list(like="^A3")))
  expect_equal(nrow(tt),  3)
})
  
context("Signatures API")
test_that("A signature can be retrieved", {  
  skip_if_devel()
  tt <- sl$clue("sigs", where_clause=list(pert_desc="sirolimus"), 
                fields = c("is_gold", "distil_id"), 
                unpack_sigs = FALSE)
  expect_equal(nrow(tt),  202)
})
test_that("Instances can be unpacked from signature endpoint", {  
  skip_if_devel()
  tt <- sl$clue("sigs", where_clause=list(pert_desc = "sirolimus"), 
                fields = c("is_gold", "distil_id"), 
                unpack_sigs = TRUE)
  expect_equal(nrow(tt),  664)
})
test_that("A list of instance ids can be retrieved", {
  skip_if_devel()
  tt <- sl$clue.instances(where_clause=list("pert_desc" = "sirolimus", 
                                            "is_gold" = TRUE, 
                                            cell_id = 'MCF7'))
  expect_equal(length(tt), 62)
})
test_that("Plate ids can be extracted from dataframe", {
  skip_if_devel()
  tt <- sl$clue("sigs", where_clause=list(pert_desc = "sirolimus"), 
                fields = c("is_gold", "distil_id"), unpack_sigs = TRUE)
  tt <- sl$distil.to.plate(tt)
  expect_equal(length(tt), 664)
  expect_equal(tt[1], "ASG001_MCF7_24H_X1")
})
test_that("Signatures can be retrieved by id", {
  skip_if_devel()
  tt <- sl$clue("sigs", where_clause=list(pert_desc = "sirolimus"),
                fields = c("is_gold", "distil_id"), unpack_sigs = TRUE)
  tt <- sl$clue("sigs", ids=tt$distil_id)
  tt <- sl$distil.to.plate(tt)
  expect_equal(length(tt), 664)
})

context("Profiles API")
test_that("Vehicle control can be retrieved for perturbations", {  
  skip_if_devel()
  tt <- sl$clue.vehicle(c("AML001_CD34_24H_X1_F1B10:A03", 
                          "ASG001_MCF7_24H_X1_B7_DUO52HI53LO:F13"))
  expect_equal(nrow(tt),  2)
  expect_equal(tt[1,3], "DMSO")
})
test_that("Instance ids can be retrieved", {
  skip_if_devel()
  tt <- sl$clue.instances(where_clause=list(pert_desc = "sirolimus"))
  expect_equal(length(tt), 664)
})
test_that("Profiles can be retrieved by id", {
  skip_if_devel()
  tt <- sl$clue("profiles", where_clause=list(pert_desc = "sirolimus"), 
                fields = c("pert_dose", "distil_id"))
  expect_equal(ncol(tt), 2)
  tt <- sl$clue("profiles", ids=tt$distil_id)
  expect_equal(length(tt), 39)
})


context("Eset")
test_that("Eset can be created", {  
  skip_if_devel()
  sl$close()
  sl <- Slinky$new(user_key, system.file("extdata", 
                                         "demo.gctx", 
                                         package="slinky"))
  eset <- sl$toEset(index = list(1:978, 1:10), 
                    info_file = system.file("extdata", 
                                            "demo_inst_info.txt", 
                                            package="slinky"))
  expect_equal(pData(eset)[1,1],  "CPC020_A375_6H_X1_B4_DUO52HI53LO:P17")
  expect_true(all.equal(exprs(eset)[5,10],  7.3669, tol=0.00001))
  expect_equal(as.numeric(nrow(eset)),  978)
  expect_equal(as.numeric(ncol(eset)),  10)
  sl$close()
})
test_that("Eset can be created by where clause", {  
  skip_if_devel()
  sl$close()
  sl <- Slinky$new(user_key,
                   gctx = system.file("extdata", "demo.gctx", package="slinky"),
                   info = system.file("extdata", "demo_inst_info.txt", package="slinky"))
  where_clause = list("sig_id"='CPC004_A375_6H:BRD-K79131256-001-08-8:10')
  eset <- sl$toEset(where_clause=where_clause)
  expect_equal(pData(eset)[1,1],  "CPC004_A375_6H_X1_B3_DUO52HI53LO:H04")
  expect_equal(as.numeric(nrow(eset)),  978)
  expect_equal(as.numeric(ncol(eset)),  3)
  sl$close()
})

## fixme
test_that("Eset can be created with controls id'd automatically", {  
  skip_if_devel()
  sl$close()
  sl <- Slinky$new(user_key,
                   gctx = system.file("extdata", "demo.gctx", package="slinky"),
                   info = system.file("extdata", "demo_inst_info.txt", package="slinky"))
  where_clause = list("pert_iname"="amoxicillin","cell_id"="A375")
  eset <- sl$toEset(where_clause=where_clause, controls=TRUE)
  expect_equal(pData(eset)[1,1],  "CPC020_A375_6H_X1_B4_DUO52HI53LO:P17")
  expect_equal(as.numeric(nrow(eset)),  978)
  expect_equal(as.numeric(ncol(eset)),  54)
  sl$close()
})

test_that("Failed query exits somewhat gracefully", {  
  skip_if_devel()
  sl$close()
  sl <- Slinky$new(user_key,
                   system.file("extdata", "demo.gctx", package="slinky"),
                   info = system.file("extdata", "demo_inst_info.txt", package="slinky"))
  expect_warning(tt <- sl$toEset(where_clause = list(pert_iname="foobar")), "no results")
  expect_null(tt)
  sl$close()
})



context("GCTX file I/O")
gctx <- system.file("extdata", "demo.gctx", package="slinky")
sl <- Slinky$new(user_key, gctx)
test_that("Column names can be retrieved", {  
  skip_if_devel()
  tt <- sl$colnames(index=list(1:3))
  expect_equal(length(tt),  3)
  expect_equal(tt[1], "CPC020_A375_6H_X1_B4_DUO52HI53LO:P17")
})
test_that("Row names can be retrieved", {  
  skip_if_devel()
  tt <- sl$rownames(index=list(1:3))
  expect_equal(length(tt),  3)
  expect_equal(tt[1], "5720")
})
test_that("Matrix data can be read", {  
  skip_if_devel()
  tt <- sl$readGCTX(index=list(1:3, 5:10))
  expect_equal(nrow(tt), 3)
  expect_equal(ncol(tt), 6)
  expect_equal(colnames(tt)[6], "CPC020_A375_6H_X1_B4_DUO52HI53LO:F05")
})
sl$close()
test_that("Matrix data can be read with file specified by constructor", {  
  skip_if_devel()
  sl = Slinky$new(key = user_key, gctx = system.file("extdata", 
                                                     "demo.gctx", 
                                                     package="slinky"))
  tt <- sl$readGCTX(gctx, index=list(1:3, 5:10))
  expect_equal(nrow(tt), 3)
  expect_equal(ncol(tt), 6)
  expect_equal(colnames(tt)[6], "CPC020_A375_6H_X1_B4_DUO52HI53LO:F05")
})

context("Characteristic direction")
sl <- Slinky$new(user_key, system.file("extdata", 
                                       "demo.gctx", 
                                       package="slinky"))
eset <- sl$toEset(index = list(1:978, 1:10), 
                  info_file = system.file("extdata", 
                                          "demo_inst_info.txt", 
                                          package="slinky"))
test_that("Characteristic direction based on expression sets", {  
  skip_if_devel()
  tt <- sl$chDir(eset[, 1:4], eset[, 5:10])
  expect_equivalent(tt[1], -0.00637839889854218029807)
})
test_that("Characteristic direction can be calculated based on two matrices", {  
  skip_if_devel()
  tt <- sl$chDir(exprs(eset[, 1:4]), exprs(eset[, 5:10]))
  expect_equivalent(tt[1], -0.00637839889854218029807)
})

