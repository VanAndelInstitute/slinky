
#' @export
#' @importFrom foreach %dopar%
#' @importFrom dplyr %>% mutate
#' @importFrom tidyr unnest
Slinky$methods(clue = function(endpoint=c("sigs", "cells", "genes", "perts", "plates", "profiles", "rep_drugs", "pcls"), 
                               fields="", 
                               where_clause="", 
                               limit=0, 
                               count=FALSE, 
                               unpack_sigs = TRUE,
                               cl=NULL) {
  "Wrapper for clue.io REST calls
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{endpoint} The endpoint to query.  Default is 'sigs'.}
  \\item{\\code{fields} Optional vector of fields to return.}
  \\item{\\code{where_clause} Optional where_clause clause.  Must be named list (e.g. list(field='value')}
  \\item{\\code{unpack_sigs} The sigs endpoint returns multiple distil_ids per row.  Should we unpack these to one per row?}
  \\item{\\code{cl} Optional cluster object to parallelize data retrieval.}
  }}
  \\subsection{Return Value}{Data returned by Slinky.api as a data.frame}"
  endpoint = match.arg(endpoint)
  
  key <- .self$.user_key
  base <- .self$.base
  if(count) {
    return(.self$clue.count(endpoint, where_clause))
  } else {
    count <- .self$clue.count(endpoint, where_clause)
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
  
  if(endpoint == "sigs" && unpack_sigs) {
    dat <- dat %>%
      mutate(distil_id = strsplit(distil_id, "\\|")) %>%
      unnest(distil_id)
  }
  
  dat
})

#' @export
Slinky$methods(clue.instances = function(where_clause = NULL) {
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
  tt <- .self$clue("sigs", fields = c("distil_id"), where_clause = where_clause, unpack_sigs = TRUE)
  return(as.character(unlist(sapply(tt[,1], function(x) { strsplit(x, "\\|") }))))
})

#' @export
Slinky$methods(clue.count = function(endpoint=c("sigs", "cells", "genes", "perts", "plates", "profiles", "rep_drugs", "pcls"), 
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

