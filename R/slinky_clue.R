# sanity check
Slinky$methods(
    .check.clue = function() {
        res <- httr::GET(url = "clue.io/api")
        grepl("API Access", httr::content(res, as = "text"))
    }
)

#' @export
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom dplyr %>% mutate
#' @importFrom tidyr unnest
#' @importFrom httr GET
Slinky$methods(clue = function(endpoint = c("sigs",
                                            "cells",
                                            "genes",
                                            "perts",
                                            "plates",
                                            "profiles",
                                            "rep_drugs",
                                            "pcls"),
                                fields = "",
                                where_clause = NULL,
                                ids = NULL,
                                limit = 0,
                                count = FALSE,
                                unpack_sigs = TRUE,
                                cl = NULL,
                                verbose = FALSE) {
    "Wrapper for clue.io REST calls
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{endpoint} The endpoint to query.  Default is 'sigs'.}
        \\item{\\code{fields} Optional vector of fields to return.}
        \\item{\\code{where_clause} Optional where_clause clause.  Must be
            named list (e.g. list(field='value')}
        \\item{\\code{ids} Optional vector of ids to fetch for sigs or profiles
            endpoints.  Should not be used together with where_clause or
            count.}
        \\item{\\code{unpack_sigs} The sigs endpoint returns multiple
            distil_ids per row.  Should we unpack these to one per row?}
        \\item{\\code{cl} Optional cluster object to parallelize this
            operation. If verbose is TRUE, use this pattern in order for
            progress bar to update:
            \\code{cl <- parallel::makeCluster(4, outfile=\"\")}}
        \\item{\\code{verbose} Do you want to know how things are going?
            Default is false.}
    }}
    \\subsection{Return Value}{Data returned by Slinky.api as a data.frame}"

    if (!.self$.check.clue()) {
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
    endpoint = match.arg(endpoint)

    key <- .self$.user_key
    base <- .self$.base
    if (count) {
        return(.self$clue.count(endpoint, where_clause))
    } else {
        if (length(ids)) {
            count <- length(ids)
        } else {
            count <- .self$clue.count(endpoint, where_clause)
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

    dat <- foreach::foreach(
        i = 1:ceiling(count / 1000),
        .combine = dplyr::bind_rows,
        .export = c(".self"),
        .packages = c("httr", "jsonlite")) %dopar% {

            if (verbose) {
                setTxtProgressBar(pb, i)
            }

            if (i * 1000 > limit) {
                lim_filter <- limit %% 1000
            }

            ii <- (i - 1) * 1000

            if (length(ids)) {
                where_clause = list(distil_id =
                                    list(inq = c(ids[ii:(ii + lim_filter)])))
                skip <- 0
            } else {
                skip <- ii
            }

            res <- httr::GET(
                url = .self$.base,
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
        dat <- dat %>%
            mutate(distil_id = strsplit(distil_id, "\\|")) %>%
            unnest(distil_id)
    }
    dat
})


#' @export
Slinky$methods(clue.vehicle = function(ids, verbose = FALSE) {
    "Fetch the vehicle control applicable to given ids (distil_id).  Expects
        that perturbagen is of type trt_cp.
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{ids} The distil_id(s) to lookup.}
        \\item{\\code{verbose} Do you want to know how things are going?
            Default is FALSE.}
    }}
    \\subsection{Return Value}{The name of the vehicle control for the queried
        perturbagen(s).}
    \\subsection{Details}{This is a convenience wrapper to the profiles API
        which queries clue.io and unwraps response. }"

    if (length(ids) > 1) {
        where_clause = list(distil_id = list(inq = c(ids)))
    } else {
        where_clause = list(distil_id = ids)
    }
    res <- .self$clue(
        "profiles",
        fields = c("pert_id",
                    "pert_iname", "pert_vehicle"),
        where_clause = where_clause,
        verbose = verbose
    )
    return(res)
})



#' @export
Slinky$methods(clue.instances = function(where_clause = NULL,
                                            verbose = FALSE,
                                            poscon = c("keep", "omit")) {
    "Convenience wrapper to query function to retrieve instance ids meeting
    specified criteria.
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{where_clause} Filter terms, as a list of terms, e.g.
            \\code{list(pert_type='trt_cp', 'is_gold'=TRUE)}.  Terms will be
            joined by AND logic.}
        \\item{\\code{verbose} Do you want to know how things are going?
            Default is false.}
        \\item{\\code{poscon} Instances of type \\code{trt_poscon} are recoded
            as \\code{trt_cp} in clue.io's \\code{sigs} endpoint.
            This can lead to unexpected results downstream.  To drop these
            instances, specify \\code{poscon='omit'}}
    }}
    \\subsection{Return Value}{Vector of ids matching criteria.}
    \\subsection{Details}{This is a convenience wrapper to the signature API
        which queries clue.io and unwraps response. }"

    if (!length(where_clause) && !length(ids)) {
        stop("The 'where_clause' argument must be specified.")
    }
    poscon = match.arg(poscon)
    ids <- .self$clue(
        "sigs",
        fields = c("distil_id"),
        where_clause = where_clause,
        unpack_sigs = TRUE,
        verbose = verbose
    )
    ids <- as.character(unlist(sapply(ids[, 1], function(x) {
        strsplit(x, "\\|")
    })))

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

#' @export
Slinky$methods(clue.count = function(endpoint = c("sigs",
                                                    "cells",
                                                    "genes",
                                                    "perts",
                                                    "plates",
                                                    "profiles",
                                                    "rep_drugs",
                                                    "pcls"),
                                    where_clause = "") {
    "Wrapper for Slinky.io REST calls to retrieve record counts
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{endpoint} The endpoint to query, default is 'sigs'.}
        \\item{\\code{where_clause} Optional where_clause clause.  Must be
        named list (e.g. list(field='value')}
    }}
    \\subsection{Return Value}{Count of records satisfying query}"

    endpoint = match.arg(endpoint)
    key <- .self$.user_key
    base <- .self$.base
    if (class(where_clause) == "list") {
        query = list(
            where = jsonlite::toJSON(where_clause,
                                    auto_unbox = TRUE),
            user_key = key
        )
    } else {
        query = list(user_key = key)
    }
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

