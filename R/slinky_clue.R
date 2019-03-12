# sanity check
.checkClue = function() {
  res <- httr::GET(url = "clue.io/api")
  grepl("API Access", httr::content(res, as = "text"))
}

#' clue
#' 
#' Wrapper for clue.io REST calls
#' 
#' @param x a Slinky Object
#' @param endpoint  The endpoint to query.  Default is 'sigs'.
#' @param fields  Optional vector of fields to return.
#' @param where_clause  Optional where_clause clause.  Must be
#'    named list (e.g. list(field='value')
#' @param ids  Optional vector of ids to fetch for sigs or profiles
#'    endpoints.  Should not be used together with where_clause or
#'    count.
#' @param unpack_sigs  The sigs endpoint returns multiple
#'    distil_ids per row.  Should we unpack these to one per row?
#' @param poscon  Instances of type \code{trt_poscon} are recoded
#'    as \code{trt_cp} in clue.io's \code{sigs} endpoint.
#'    This can lead to unexpected results downstream.  To keep these
#'    instances, specify \code{poscon='keep'} 
#' @param limit  Optional limit to number of instances (samples to return
#' @param count  Should we just return the count of intances satisfying 
#'    the query rather than the data? Default is FALSE.
#' @param cl  Optional cluster object to parallelize this
#'    operation. If verbose is TRUE, use this pattern in order for
#'    progress bar to update:
#'    \code{cl <- parallel::makeCluster(4, outfile=\"\")}
#' @param verbose  Do you want to know how things are going?
#'    Default is false.
#' @return Data returned by Slinky.api as a data.frame
#' @name clue
#' @rdname clue
#' @import foreach
#' @importFrom dplyr %>% mutate
#' @importFrom tidyr unnest
#' @importFrom httr GET
#' @examples
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package = 'slinky'))
#' amox <- clue(sl, where_clause = list("pert_iname" = "amoxicillin", 
#'                                     "cell_id" = "MCF7",
#'                                     "is_gold" = TRUE))
setGeneric("clue",
	function(x, endpoint = c("sigs",
                                            "cells",
                                            "genes",
                                            "perts",
                                            "plates",
                                            "profiles",
                                            "rep_drugs",
	                                          "rep_drug_indications",
	                                          "pcls"),
                               fields = "",
                               where_clause = NULL,
                               ids = NULL,
                               limit = 0,
                               count = FALSE,
                               unpack_sigs = TRUE,
	                             poscon = c("omit", "keep"),
                               cl = NULL,
                               verbose = FALSE) {
	standardGeneric("clue")
	}
)
#' @rdname clue
#' @exportMethod clue
#' @importFrom dplyr bind_rows %>%
#' @aliases clue,Slinky-method
setMethod("clue", signature(x = "Slinky"),
function(x, endpoint = c("sigs",
                                            "cells",
                                            "genes",
                                            "perts",
                                            "plates",
                                            "profiles",
                                            "rep_drug_indications",
                                            "rep_drugs",
                                            "pcls"),
                               fields = "",
                               where_clause = NULL,
                               ids = NULL,
                               limit = 0,
                               count = FALSE,
                               unpack_sigs = TRUE,
                               poscon = c("omit", "keep"),
                               cl = NULL,
                               verbose = FALSE) 
{

  
  if (!.checkClue()) {
    stop("Could not connect to clue.io APi.  Please verify connecivity.")
  }
  
  if (sum(length(ids) & count & length(where_clause)) > 1) {
    stop(
      "In call to Slinky$clue, the ids, where_clause, and count arguments
      are mutually exclusive.  Ensure only 1 provided."
    )
  }
  
  if (length(ids) && !endpoint %in% c("sigs", "profiles")) {
    stop("Specifying id list currently only supported for sigs and ",
         "profiles endpoints")
  }
  endpoint <- match.arg(endpoint)
  poscon <- match.arg(poscon)
  key <- x@user_key
  base <- x@base
  
  if (count) {
    return(clueCount(x, endpoint, where_clause))
  } else {
    if (length(ids)) {
      count <- length(ids)
    } else {
      count <- clueCount(x, endpoint, where_clause)
    }
  }
  
  if (count == 0)
    return(NULL)
  
  if (limit > 0) {
    count <- min(limit, count)
  } else {
    limit <- count
  }
  
  if (length(cl)) {
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  
  # we cannot retrieve more than 1000 records at a time, nor can we send
  # more than 1000 ids in an inq statement so throttle to 1000 records per
  # request, unless user has requested fewer
  if (limit < 1000) {
    lim_filter <- limit
  } else {
    lim_filter <- 1000
  }
  if (verbose) {
    message("Retrieving metadata from clue.io")
    pb <- txtProgressBar(
      min = 0,
      max = ceiling(count / 1000),
      initial = 0,
      width = 100,
      style = 3
    )
  }
  
  if (poscon == "omit") {
    pc <- readRDS(system.file("extdata", "trt_poscon.rds",
                              package = "slinky"))
  }
  
  i <- 0 # keep R CHECK happy.
  
  dat <- foreach(
    i = seq_len(ceiling(count / 1000)),
    .combine = dplyr::bind_rows,
    .export = c("x"),
    .packages = c("httr", "jsonlite")) %dopar% {
      
      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }
      
      if (i * 1000 > limit) {
        lim_filter <- limit %% 1000
      }
      
      ii <- (i - 1) * 1000
      
      if (length(ids)) {
        where_clause = list("distil_id" =
                              list(inq = c(ids[seq(ii, ii + lim_filter)])))
        skip <- 0
      } else {
        skip <- ii
      }
      if (poscon == "omit") {
        ix <- which(ids %in% pc$distil_id)
        if (length(ix)) {
          ids <- ids[-ix]
        }
      }
      
      res <- httr::GET(
        url = x@base,
        path = paste0("api/", endpoint, "/"),
        query = list(
          filter = jsonlite::toJSON(
            list(
              fields = fields,
              limit = lim_filter,
              where = where_clause,
              skip = skip
            ),
            auto_unbox = TRUE
          ),
          user_key = key
        )
      )
      
      # flatten the data structures to 2D
      fl <- lapply(httr::content(res, as = "parsed"),
                   function(x) {
                     lapply(x, unlist)
                   })
      fl <- lapply(fl, function(x) {
        lapply(x, paste, collapse = "|")
      })
      dplyr::bind_rows(lapply(fl,
                              as.data.frame,
                              stringsAsFactors = FALSE))
    }
  if (length(dat) == 0) {
    return(NULL)
  }
  if (endpoint == "sigs" && unpack_sigs) {
    distil_id <- NULL # prevent "no visible binding" on R CMD CHECK
    dat <- dat %>%
      mutate("distil_id" = strsplit(distil_id, "\\|")) %>%
      unnest(distil_id)
  }
  dat
  })


#' clueVehicle
#' 
#' Fetch the vehicle control applicable to given ids (distil_id).  Expects
#'    that perturbagen is of type trt_cp.
#'    
#' @param x a Slinky Object
#' @param ids  The distil_id(s) to lookup.
#' @param verbose  Do you want to know how things are going?
#'    Default is FALSE.
#' @return The name of the vehicle control for the queried
#'    perturbagen(s).
#'    This is a convenience wrapper to the profiles API
#'    which queries clue.io and unwraps response. 
#' @name clueVehicle
#' @rdname clueVehicle
#' @examples
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package = 'slinky'))
#' amox <- clue(sl, where_clause = list("pert_iname" = "amoxicillin", 
#'                                     "cell_id" = "MCF7",
#'                                     "is_gold" = TRUE))
#' amox.ctrl <- clueVehicle(sl, amox$distil_id)
#' 
setGeneric("clueVehicle",
	function(x, ids, verbose = FALSE) {
	standardGeneric("clueVehicle")
	}
)
#' @rdname clueVehicle
#' @exportMethod clueVehicle
#' @aliases clueVehicle
setMethod("clueVehicle", signature(x = "Slinky"),
function(x, ids, verbose = FALSE) 
{

  
  if (length(ids) > 1) {
    where_clause = list(distil_id = list(inq = c(ids)))
  } else {
    where_clause = list(distil_id = ids)
  }
  res <- clue(x, 
    "profiles",
    fields = c("pert_id",
               "pert_iname", "pert_vehicle"),
    where_clause = where_clause,
    verbose = verbose
  )
  return(res)
})



#' clueInstances
#' 
#' Convenience wrapper to query function to retrieve instance ids meeting
#'    specified criteria.
#'    
#' @param x a Slinky Object
#' @param where_clause  Filter terms, as a list of terms, e.g.
#'    \code{list(pert_type='trt_cp', 'is_gold'=TRUE)}.  Terms will be
#'    joined by AND logic.
#' @param verbose  Do you want to know how things are going?
#'    Default is false.
#' @param poscon  Instances of type \code{trt_poscon} are recoded
#'    as \code{trt_cp} in clue.io's \code{sigs} endpoint.
#'    This can lead to unexpected results downstream.  To keep these
#'    instances, specify \code{poscon='keep'} 
#' @return Vector of ids matching criteria.
#'    This is a convenience wrapper to the signature API
#'    which queries clue.io and unwraps response. 
#' @name clueInstances
#' @rdname clueInstances
setGeneric("clueInstances",
	function(x, where_clause = NULL,
                                        verbose = FALSE,
                                        poscon = c("omit", "keep")) {
	standardGeneric("clueInstances")
	}
)
#' @examples
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package = 'slinky'))
#' amox_ids <- clueInstances(sl, where_clause = list("pert_iname" = "amoxicillin", 
#'                                     "cell_id" = "MCF7",
#'                                     "is_gold" = TRUE))
#' 
#' @rdname clueInstances
#' @exportMethod clueInstances
#' @aliases clueInstances
setMethod("clueInstances", signature(x = "Slinky"),
function(x, where_clause = NULL,
                                        verbose = FALSE,
                                        poscon = c("omit", "keep")) 
{

  
  if (!length(where_clause) && !length(ids)) {
    stop("The 'where_clause' argument must be specified.")
  }
  poscon = match.arg(poscon)
  ids <- clue(x,
    "sigs",
    fields = c("distil_id"),
    where_clause = where_clause,
    unpack_sigs = TRUE,
    verbose = verbose
  )
  ids <- as.character(unlist(vapply(ids[, 1], function(x) {
    strsplit(x, "\\|")
  }, list("a"))))
  
  if (poscon == "omit") {
    pc <- readRDS(system.file("extdata", "trt_poscon.rds",
                              package = "slinky"))
    ix <- which(ids %in% pc$distil_id)
    if (length(ix)) {
      ids <- ids[-ix]
    }
  }
  ids
})

#' clueCount
#' 
#' Wrapper for Slinky.io REST calls to retrieve record counts
#' 
#' @param x a Slinky Object
#' @param endpoint  The endpoint to query, default is 'sigs'.
#' @param where_clause  Optional where_clause clause.  Must be
#'    named list (e.g. list(field='value')
#' @return Count of records satisfying query
#' @name clueCount
#' @rdname clueCount
setGeneric("clueCount",
	function(x, endpoint = c("sigs",
                                                 "cells",
                                                 "genes",
                                                 "perts",
                                                 "plates",
                                                 "profiles",
                                                 "rep_drugs",
                      	                         "rep_drug_indications",
                                                 "pcls"),
                                    where_clause = "") {
	standardGeneric("clueCount")
	}
)
#' @examples
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package = 'slinky'))
#' amox_count <- clueCount(sl, where_clause = list("pert_iname" = "amoxicillin", 
#'                                     "cell_id" = "MCF7",
#'                                     "is_gold" = TRUE))
#'                                     
#' @rdname clueCount
#' @exportMethod clueCount
#' @aliases clueCount
setMethod("clueCount", signature(x = "Slinky"),
function(x, endpoint = c("sigs",
                                                 "cells",
                                                 "genes",
                                                 "perts",
                                                 "plates",
                                                 "profiles",
                                                 "rep_drugs",
                                                 "rep_drug_indications",
                                                 "pcls"),
                                    where_clause = "") 
{

  
  endpoint = match.arg(endpoint)
  key <- x@user_key
  base <- x@base
  if (is(where_clause, "list")) {
    query = list(
      where = jsonlite::toJSON(where_clause,
                               auto_unbox = TRUE),
      user_key = key
    )
  } else {
    query = list(user_key = key)
  }
  

  #  https://api.clue.io/api/rep_drugs/count?where=%7B%22pert_iname%22%3A%22tacrolimus%22%7D&user_key=MYSUSERKEY
  res <- httr::GET(
    url = base,
    path = paste0("api/", endpoint, "/count"),
    query = query
  )
  res <- httr::content(res)
  if (length(res$error)) {
    stop("Call to clue.io API returned this error: ", res$error)
  }
  count <- res$count
  if (!length(count)) {
    count <- 0
  }
  return(count)
})

