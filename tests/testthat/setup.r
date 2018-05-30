library(slinky)
library(testthat)


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
                          as = "parsed")$user_key
gctx <- system.file("extdata", "demo.gctx", package = "slinky")
info = system.file("extdata",
                   "demo_inst_info.txt",
                   package = "slinky")
sl <- Slinky(user_key, gctx, info)

