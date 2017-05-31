#' A Reference Class to encapsulate LINCS methods
#'
#' @field .user_key private field populated at instantiation
#' @field .base private field populated at instantiation
Slinky <- setRefClass("Slinky",
                      fields = list(.user_key = "character",
                                    .base = "character"),
                      methods = list(
                      )
)

Slinky$methods(initialize = function(user_key=NULL) {

  if(!length(user_key)) {
    user_key = Sys.getenv("CLUE_API_KEY")
    if(nchar(user_key) < 10) {
      stop("No user key provided.  Either provide user_key argument or set Slinky_API_KEY env variable")
    }
  }
  .self$.user_key = user_key
  .self$.base = "https://api.clue.io"
})

Slinky$methods(count = function(path, where_clause="") {
  "Wrapper for Slinky.io REST calls to retrieve record counts
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{path} The endpoint to query.}
  \\item{\\code{fields} Optional vector of fields to return.}
  \\item{\\code{where_clause} Optional where_clause clause.  Must be named list (e.g. list(field='value')}
  }}
  \\subsection{Return Value}{Count of records satisfying query}"

  user_key <- .self$.user_key
  base <- .self$.base
  if(class(where_clause) == "list") {
    query = list(where=jsonlite::toJSON(where_clause, auto_unbox = TRUE), user_key=user_key)
  } else{
    query = list(user_key=user_key)
  }
  res <- GET(url = base, 
             path=paste0("api/", path, "/count"), 
             query = query)
  count <- content(res)$count
  if(!length(count)) {
    count <- 0
  }
  return(count)
})

Slinky$methods(fetch = function(path, fields="", where_clause="", limit=0, count=FALSE, cl=NULL) {
  "Wrapper for Slinky.io REST calls
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{path} The endpoint to query.}
  \\item{\\code{fields} Optional vector of fields to return.}
  \\item{\\code{where_clause} Optional where_clause clause.  Must be named list (e.g. list(field='value')}
  \\item{\\code{cl} Optional cluster object to parallelize data retrieval.}
  }}
  \\subsection{Return Value}{Data returned by Slinky.api as a data.frame}"
  user_key <- .self$.user_key
  base <- .self$.base
  if(count) {
    return(.self$count(path, where_clause))
  } else {
    count <- .self$count(path, where_clause)
  }
  
  if(count==0) return(NULL)

  if(limit > 0) {
    count <- min(limit, count)
  }
  
  if(length(cl)) {
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }

  if(limit > 0 & limit < 1000) {
    lim_filter <- limit
  } else {
    lim_filter <- 1000
  }

  dat <- foreach(i=1:ceiling(count/1000), 
                 .combine = rbind.fill,
                 .export=c(".self"),
                 .packages = c("httr", "jsonlite")) %dopar% 
  {
    if(i * 1000 > limit & limit > 0) {
      lim_filter <- limit %% 1000
    }
    res <- GET(url = .self$.base, 
             path=paste0("api/", path, "/"), 
             query = list(filter=jsonlite::toJSON(list(fields=fields, 
                                                       limit=lim_filter,
                                                       where=where_clause, 
                                                       skip=(i-1)*1000), 
                                                  auto_unbox = TRUE), 
                          user_key=user_key))
    # flatten the data structures to 2D
    fl <- lapply(content(res, as = "parsed"), function(x) { lapply(x, unlist) })
    fl <- lapply(fl, function(x) { lapply(x, paste, collapse="|")})
    rbind.fill(lapply(fl, as.data.frame, stringsAsFactors=FALSE))
  }
  if(limit > 0 & nrow(dat) > limit) {
    dat<-dat[1:limit,]
  }
  dat
})

Slinky$methods(gctx.colnames = function(file, index=NULL) {
  "Retrieve column names from LINCS gctx datafile
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Path to gctx file.}
  \\item{\\code{index} Optional list providing which colnames to return.}
  }}
  \\subsection{Return Value}{Names of columns from gctx file}
  \\subsection{Details}{The gctx file is an HDF5 formatted file with several sections (groups) containing the column and
  row level metadata as well as the expression data itself.  If `index` is provided, it should be a list of extent one providing 
  a vector of column indices for which to retrn ids.  For example, index=list(c(1,2,3,10)).}"
  if((length(index) & !typeof(index) == "list") || length(index) > 1) {
    stop("Index must be a list of extent 1 provided indices of instanceIds to return")
  }
  if(!length(index)) {
    info <- h5dump(file, load=FALSE)
    index <- list(1:info$`0`$META$COL$id$dim)    
  }
  h5read(file, name = "0/META/COL/id", index=index)
})

Slinky$methods(gctx.rownames = function(file, index=NULL) {
  "Retrieve rown names from LINCS gctx datafile
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Path to gctx file.}
  \\item{\\code{index} Optional list providing which rownames to return.}
  }}
  \\subsection{Return Value}{Names of rows from gctx file}
  \\subsection{Details}{The gctx file is an HDF5 formatted file with several sections (groups) containing the column and
  row level metadata as well as the expression data itself.  If `index` is provided, it should be a list of extent one providing 
  a vector of row indices for which to retrn ids.  For example, index=list(c(1,2,3,10)).}"
  if((length(index) & !typeof(index) == "list") || length(index) > 1) {
    stop("Index must be a list of extent 1 provided indices of instanceIds to return")
  }
  if(!length(index)) {
    info <- h5dump(file, load=FALSE)
    index <- list(1:info$`0`$META$ROW$id$dim)    
  }
  h5read(file, name = "0/META/ROW/id", index=index)
})

Slinky$methods(gctx.read = function(file, index=NULL) {
  "Read portions of data matrix from LINCS gctx datafile
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Path to gctx file.}
  \\item{\\code{index} list of extent 2 providing which columns and rows to return.  Since the LINCS datasets are relatively
  massive, the index argument is required to avoid inadvertently trying to slurp the entire dataset into memory.}
  }}
  \\subsection{Return Value}{Matrix of expression data with rownames and colnames appropriately set}"
  if(!length(index)) {
    stop("You must provide index argument.  E.g. index=list(1:10, 1:5)")
  }
  data <- h5read(file, name = "0/DATA/0/matrix", index=index)
  colnames(data) <- .self$gctx.colnames(file, index=list(index[[2]]))
  rownames(data) <- .self$gctx.rownames(file, index=list(index[[1]]))
  data
})

Slinky$methods(score = function(data, method, ...) {
  "Provide a uniform interface to scoring functions.
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{data} Matrix of expression data to score. }
  \\item{\\code{method} Scoring method to use.  Only \\code{ks} is presently supported.}
  \\item{\\code{...} Additional arguments for \\code{method}.}
  }}
  \\subsection{Return Value}{Vectors of scores, one per column of data}"
  if(method == "ks") {
    return(.self$.ks(data, ...))
  } else if(method == "xsum") {
    return(.self$.xsum(data, ...))
  } else {
    return(NULL)
  }
})

Slinky$methods(.ks = function(data, up, down, plot=FALSE) {
  "Calculate KS based enrichment statistic (CMAP method, Lamb, et al.)
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{data} Matrix of expression data to score. }
  \\item{\\code{up} List of upregulated gene ids (must match rownames(data)).}
  \\item{\\code{down} List of downregulated gene ids (must match rownames(data)).}
  \\item{\\code{plot} Whether to plot GSEA style plot for each column (not implemented).}
  }}
  \\subsection{Return Value}{Vectors of scores, all in (0,1), one per column of data}"
  f <- function(ix, n) {
    sc <- rep(-1/(n-length(ix)), n)
    # account for ties
    ix <- round(ix)
    sc[ix] <- 0
    for(i in ix) {
      sc[i] <- sc[i] + 1/length(ix) 
    }
    sc <- cumsum(sc) 
    if(max(sc) > -min(sc)) {
      return(max(sc))
    } else {
      return(min(sc))
    }
  }
  data <- apply(-data, 2, rank)
  up <- which(rownames(data) %in% up)
  down <- which(rownames(data) %in% down)
  apply(data, 2, function(x) {
    u <- f(x[up], length(x))
    d <- f(x[down], length(x))
    if(sign(u) == sign(d)) {
      return(0)
    } else {
      return((0.5*(abs(u) + abs(d))) * sign(u))
    }
  })
})

Slinky$methods(.xsum = function(data, up, down, n=100) {
  "Calculate Extreme Sum stat (Agerwal et al.))
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{data} Matrix of expression data to score. }
  \\item{\\code{up} List of upregulated gene ids (must match rownames(data)).}
  \\item{\\code{down} List of downregulated gene ids (must match rownames(data)).}
  \\item{\\code{n} How many genes to consider for scoring among up and down regulated genes in each column of the data matrix.}
  }}
  \\subsection{Return Value}{Vectors of scores, all in one per column of data}"
  
  up.ix <- which(rownames(data) %in% up)
  down.ix <- which(rownames(data) %in% down)
  f <- function(a) {
    a_r <- rank(a)
    changed <- a * (a_r > ( length(a_r) - n) | a_r < n)
    sum(changed[up.ix]) - sum(changed[down.ix], na.rm=TRUE)  
  }
  apply(data, 2, f)
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
