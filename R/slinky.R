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
                 .packages = c("httr", "jsonlite")) %do% 
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
