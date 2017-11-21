#' A Reference Class to encapsulate LINCS methods
#'
#' @field .user_key private field populated at instantiation
#' @field .base private field populated at instantiation
#' @field .gctx private field optionally populated at instantiation specifying path to gctx file
#' @import methods
#' @export Slinky
#' @exportClass Slinky
Slinky <- setRefClass("Slinky",
                      fields = list(.user_key = "character",
                                    .base = "character",
                                    .gctx = "ANY"), # allow null
                     methods = list(
                                     )
)

#' @export
Slinky$methods(initialize = function(key=NULL, gctx=NULL) {
  "Create a Slinky object
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{key} Your clue.io user key.  If not specified, CLUE_API_KEY environment variable must be set.}
  \\item{\\code{gctx} Optional path to gctx file containing data you want to work with.  Can be specified later if desired.}
  }}
  \\subsection{Return Value}{Slinky object}"
  
  if(!length(key)) {
    key = Sys.getenv("CLUE_API_KEY")
    if(nchar(key) < 10) {
      stop("No user key provided.  Either provide user_key argument or set CLUE_API_KEY env variable")
    }
  }
  .self$.user_key = key
  .self$.base = "https://api.clue.io"
  if(length(gctx)) {
    .self$.gctx = gctx;
  } else {
    .self$.gctx = NULL;
  }
})

#' @export
Slinky$methods(count = function(endpoint=c("sigs", "cells", "genes", "perts", "plates", "profiles", "rep_drugs", "pcls"), 
                                where_clause="") {
  "Wrapper for Slinky.io REST calls to retrieve record counts
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{endpoint} The endpoint to query, default is 'sigs'.}
  \\item{\\code{where_clause} Optional where_clause clause.  Must be named list (e.g. list(field='value')}
  }}
  \\subsection{Return Value}{Count of records satisfying query}"

  endpoint = match.arg(endpoint)
  key <- .self$.user_key
  base <- .self$.base
  if(class(where_clause) == "list") {
    query = list(where=jsonlite::toJSON(where_clause, auto_unbox = TRUE), user_key=key)
  } else{
    query = list(user_key=key)
  }
  res <- httr::GET(url = base, 
             path=paste0("api/", endpoint, "/count"), 
             query = query)
  count <- httr::content(res)$count
  if(!length(count)) {
    count <- 0
  }
  return(count)
})

#' @export
#' @importFrom foreach %dopar%
Slinky$methods(query = function(endpoint=c("sigs", "cells", "genes", "perts", "plates", "profiles", "rep_drugs", "pcls"), 
                                fields="", 
                                where_clause="", 
                                limit=0, 
                                count=FALSE, 
                                cl=NULL) {
  "Wrapper for Slinky.io REST calls
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{endpoint} The endpoint to query.  Default is 'sigs'.}
  \\item{\\code{fields} Optional vector of fields to return.}
  \\item{\\code{where_clause} Optional where_clause clause.  Must be named list (e.g. list(field='value')}
  \\item{\\code{cl} Optional cluster object to parallelize data retrieval.}
  }}
  \\subsection{Return Value}{Data returned by Slinky.api as a data.frame}"
  endpoint = match.arg(endpoint)
  
  key <- .self$.user_key
  base <- .self$.base
  if(count) {
    return(.self$count(endpoint, where_clause))
  } else {
    count <- .self$count(endpoint, where_clause)
  }
  
  if(count==0) return(NULL)

  if(limit > 0) {
    count <- min(limit, count)
  }
  
  if(length(cl)) {
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  if(limit > 0 & limit < 1000) {
    lim_filter <- limit
  } else {
    lim_filter <- 1000
  }

  dat <- foreach::foreach(i=1:ceiling(count/1000), 
                 .combine = plyr::rbind.fill,
                 .export=c(".self"),
                 .packages = c("httr", "jsonlite")) %dopar% 
  {
    if(i * 1000 > limit & limit > 0) {
      lim_filter <- limit %% 1000
    }
    res <- httr::GET(url = .self$.base, 
             path=paste0("api/", endpoint, "/"), 
             query = list(filter=jsonlite::toJSON(list(fields=fields, 
                                                       limit=lim_filter,
                                                       where=where_clause, 
                                                       skip=(i-1)*1000), 
                                                  auto_unbox = TRUE), 
                          user_key=key))
    # flatten the data structures to 2D
    fl <- lapply(httr::content(res, as = "parsed"), function(x) { lapply(x, unlist) })
    fl <- lapply(fl, function(x) { lapply(x, paste, collapse="|")})
    plyr::rbind.fill(lapply(fl, as.data.frame, stringsAsFactors=FALSE))
  }
  if(limit > 0 & nrow(dat) > limit) {
    dat<-dat[1:limit,]
  }
  dat
})

#' @export
Slinky$methods(query.instances = function(where_clause = NULL) {
  "Convenience wrapper to query function to retrieve instance ids meeting specified criteria
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{where_clause} Filter terms, as a list of terms, e.g. list(pert_type='trt_cp', 'is_gold'=TRUE).  Terms will be joined by AND logic.}
  }}
  \\subsection{Return Value}{Vector of ids matching criteria.}
  \\subsection{Details}{This is a convenience wrapper to the signature API which queries clue.io and unwraps response. }"
  if(!length(where_clause)) {
    stop("The 'where' argument must be specified.")
  }
  tt <- .self$query("sigs", fields = c("distil_id"), where_clause=where_clause)
  return(as.character(unlist(sapply(tt[,1], function(x) { strsplit(x, "\\|") }))))
})

#' @export
Slinky$methods(colnames = function(file=NULL, index=NULL) {
  "Retrieve column names from LINCS gctx datafile
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Path to gctx file.  May be omitted if specified when Slinky object is created.}
  \\item{\\code{index} Optional list providing which colnames to return.}
  }}
  \\subsection{Return Value}{Names of columns from gctx file}
  \\subsection{Details}{The gctx file is an HDF5 formatted file with several sections (groups) containing the column and
  row level metadata as well as the expression data itself.  If `index` is provided, it should be a list of extent one providing 
  a vector of column indices for which to retrn ids.  For example, index=list(c(1,2,3,10)).}"
  if((length(index) & !typeof(index) == "list") || length(index) > 1) {
    stop("Index must be a list of extent 1 provided indices of instanceIds to return")
  }
  if(!length(file) && !length(.self$.gctx)) {
    stop("You must specify path to gctx file.")
  }
  if(!length(file)) file = .self$.gctx
  if(!length(index)) {
    info <- rhdf5::h5dump(file, load=FALSE)
    index <- list(1:info$`0`$META$COL$id$dim)    
  }
  rhdf5::h5read(file, name = "0/META/COL/id", index=index)
})

#' @export
Slinky$methods(rownames = function(file=NULL, index=NULL) {
  "Retrieve rown names from LINCS gctx datafile
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Path to gctx file.  May be omitted if specified when Slinky object is created.}
  \\item{\\code{index} Optional list providing which rownames to return.}
  }}
  \\subsection{Return Value}{Names of rows from gctx file}
  \\subsection{Details}{The gctx file is an HDF5 formatted file with several sections (groups) containing the column and
  row level metadata as well as the expression data itself.  If `index` is provided, it should be a list of extent one providing 
  a vector of row indices for which to retrn ids.  For example, index=list(c(1,2,3,10)).}"
  if((length(index) & !typeof(index) == "list") || length(index) > 1) {
    stop("Index must be a list of extent 1 provided indices of instanceIds to return")
  }
  if(!length(file) && !length(.self$.gctx)) {
    stop("You must specify path to gctx file.")
  }
  if(!length(file)) file = .self$.gctx
  if(!length(index)) {
    info <- rhdf5::h5dump(file, load=FALSE)
    index <- list(1:info$`0`$META$ROW$id$dim)    
  }
  rhdf5::h5read(file, name = "0/META/ROW/id", index=index)
})

#' @export
Slinky$methods(readGCTX = function(file=NULL, index=NULL) {
  "Read portions of data matrix from LINCS gctx datafile
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Path to gctx file.  May be omitted if specified when Slinky object is created.}
  \\item{\\code{index} list of extent 2 providing which columns and rows to return.  Since the LINCS datasets are relatively
  massive, the index argument is required to avoid inadvertently trying to slurp the entire dataset into memory.}
  }}
  \\subsection{Return Value}{Matrix of expression data with rownames and colnames appropriately set}
  \\subsection{Details}{Failure to close connections will result in warning messages at garbage collection.  
  Calling close prior to removing object will close connections and prevent these messages.}"
  if(!length(file) && !length(.self$.gctx)) {
    stop("You must specify path to gctx file.")
  }
  if(!length(index)) {
    stop("You must provide index argument.  E.g. index=list(1:10, 1:5)")
  }
  if(!length(file)) file = .self$.gctx
  data <- rhdf5::h5read(file, name = "0/DATA/0/matrix", index=index)
  colnames(data) <- .self$colnames(file, index=list(index[[2]]))
  rownames(data) <- .self$rownames(file, index=list(index[[1]]))
  data
})

#' @export
Slinky$methods(close = function() {
  "Close any open HDF5 (gctx) file connections.
  \\subsection{Return Value}{None.  Called for side effect of closing connections.}"
  rhdf5::H5close()
})
  

# do.call(f, as.list(...))

#https://api.clue.io/api/profiles?filter={%22where%22:{%22distil_id%22:%22ASG001_MCF7_24H_X1_B7_DUO52HI53LO:A17%22}}&user_key=fe5dfc18a374f94119ea6847c439c2d2
#"pert_vehicle": "DMSO",

## Test gctx file containing tiny subset of LINCS data created as follows
#
# meta.col.id <- h5read("/media/primary/data2/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx", 
#                       name = "0/META/COL/id", 
#                       index=list(1:10))
# meta.row.id <- h5read("/media/primary/data2/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx", 
#                       name = "0/META/ROW/id", 
#                       index=list(1:5))
# data.0.matrix <- h5read("/media/primary/data2/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx", 
#                         name = "0/DATA/0/matrix", 
#                         index=list(1:5, 1:10))
# h5createFile("inst/extdata/demo.gctx")
# h5createGroup("inst/extdata/demo.gctx","0")
# h5createGroup("inst/extdata/demo.gctx","0/META")
# h5createGroup("inst/extdata/demo.gctx","0/META/COL")
# h5createGroup("inst/extdata/demo.gctx","0/META/ROW")
# h5createGroup("inst/extdata/demo.gctx","0/DATA")
# h5createGroup("inst/extdata/demo.gctx","0/DATA/0")
# 
# h5createDataset("inst/extdata/demo.gctx", "0/META/COL/id", 
#                 c(10), 
#                 size=36, 
#                 storage.mode = "character")
# h5write(meta.col.id, file="inst/extdata/demo.gctx", name="0/META/COL/id")
# 
# h5createDataset("inst/extdata/demo.gctx", "0/META/ROW/id", 
#                 c(5), 
#                 storage.mode = "integer")
# h5write(meta.row.id, file="inst/extdata/demo.gctx", name="0/META/ROW/id")
# 
# h5createDataset("inst/extdata/demo.gctx", "0/DATA/0/matrix", 
#                 c(5,10), 
#                 storage.mode = "double")
# h5write(data.0.matrix, file="inst/extdata/demo.gctx", name="0/DATA/0/matrix", 
#         file="inst/extdata/demo.gctx", 
#         name="0/META/ROW/id")
# 
