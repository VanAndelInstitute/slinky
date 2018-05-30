#' as("Slinky", "SummarizedExperiment")
#' 
#' Create SummarizedExperiment object from Slinky object.  Data will be 
#'     loaded from the GCTX file, combined with metadata from the info file, 
#'     and wrapped in a SummarizedExperiment object.  Note that this may take 
#'     a long time for the entire data set from the L1000 project.  For most 
#'     use cases, a subset will be desired (e.g. 
#'     \code{as(sl[1:987,1:50], SummarizedExperiment)}.  See the
#'     \code{\link{loadL1K}} method for a more flexible way to create a 
#'     SummarizedExperiment object based on various query parameters.
#' @name coerce
#' @rdname coerce
#' @aliases coerce,Slinky,SummarizedExperiment-method
#' @exportMethod coerce
#' @examples
#'
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package = 'slinky'))
#' sumex <- as(sl[, 1:20], "SummarizedExperiment")
#' 
#' @seealso \code{\link{loadL1K}}
#' @importFrom SummarizedExperiment SummarizedExperiment
setAs("Slinky", "SummarizedExperiment", function(from) {
  data <- readGCTX(from)
  ix <- match(base::colnames(data), from@metadata$inst_id)
  info <- from@metadata[ix, ]
  rownames(info) <- info$inst_id
  info$distil_id <- info$inst_id
  
  # this should never happen, but just a sanity check...
  if (!all.equal(base::rownames(info), 
                 as.vector(base::colnames(data)))) {
    stop("Rownames of metadata and colnames of expression data ",
         "did not match.")
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = data), 
    colData = info)
})


#' loadL1K
#' 
#' Convert data from gctx file to SummarizedExperiment, pulling metadata 
#'     from various sources.  Specifying the subset of data you want can be 
#'     done is various ways.  The simplest example is 
#'     with explicit subsetting 
#'     (e.g. \code{data <- loadL1K(sl[1:978, 1:10])}, which is equivalent to 
#'     \code{data <- as(sl[1:978, 1:10], "SummarizedExperiment")}).  However, 
#'     more sophisticated data loading can be achieved, for example by 
#'     specifing a specific perturbagen (\code{pert}), or even an explicit 
#'     \code{where_clause} which will be passed directly to the clue.io API. 
#'     Thus, this function supports users with varying degrees of familiarity 
#'     with the structure and content of the L1000 metadata.  All arguments 
#'     are optional.
#' @param x  A slinky object
#' @param ids  distil_ids to include in the Expression Set.
#' @param pert  name (pert_iname) of perturbation for which data is desired. 
#'     Supercedes value for \code{pert_iname} specified in \code{where_clause} 
#'     (if any).
#' @param cell_line  name (cell_id) of cell line for which data is desired.
#'     Supercedes value for \code{cell_id} specified in \code{where_clause} 
#'     (if any).
#' @param type  Optional type (pert_type) of perturbation for which data is desired, 
#'     should be one of \code{c("trt_cp", "trt_sh", "trt_oe")} 
#'     or NULL (default)
#' @param where_clause  Rather than specifying above terms, an explicit
#'     where_clause may be provided to identify the data to be loaded from
#'     the gctx file.  This will be passed directly the the 
#'     \href{https://clue.io/API#perts}{/pert} endpoint of 
#'     the clue.io API, and full documentation of the query options can be 
#'     reviewed at the above link.
#' @param inferred  Should the inferred (non-landmark) genes be
#'     included in the analysis? Default is TRUE. Ignored if index is
#'     specified.
#' @param gold  Should we limit to instances classified as "gold" by the L1000
#'     project (by virtue of their low inter-replicate variability)?  Default 
#'     is \code{TRUE}.
#' @param fields  Fields to include in the expression set's
#'     phenodata. Default is all available.
#' @param controls  Should same-plate controls be identified and
#'     included?  Default is FALSE.
#' @param cl  Optional cluster object to speed up data retrieval from clue.io.
#'     Please use caution...a large cluster might produce requests to the 
#'     API at an obnoxious rate.
#' @param verbose  Do you want to know how things are going?
#'     Default is FALSE
#' @return Object of type 
#'     \code{\link{SummarizedExperiment}} containing 
#'     expression and meta data.
#'     #' @name loadL1K
#' @examples
#'
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package = 'slinky'))
#' amox_gold <- clueInstances(sl, where_clause = list("pert_type" = "trt_cp",
#'                                                   "pert_iname" = "amoxicillin",
#'                                                   "cell_id" = "MCF7",
#'                                                   "is_gold" = TRUE), 
#'                           poscon = "omit")
#' ids.ctrl <- controls(sl, ids = amox_gold)$distil_id
#' amox_and_control <- loadL1K(sl, ids = c(amox_gold, 
#'                                        ids.ctrl))
#' str(SummarizedExperiment::assays(amox_and_control)[[1]]) 
#'                                       
#' @rdname loadL1K
setGeneric("loadL1K",
           function(x, 
                    ids = character(),
                    pert = character(),
                    cell_line = character(), 
                    type = NULL,
                    where_clause = character(),
                    gold = TRUE,
                    inferred = TRUE,
                    fields = character(),
                    controls = FALSE,
                    cl = NULL,
                    verbose = FALSE) {
             standardGeneric("loadL1K")
           }
)
#' @rdname loadL1K
#' @aliases loadL1K,Slinky-method
#' @exportMethod loadL1K
setMethod("loadL1K", signature(x = "Slinky"),
          function(x,
                   ids = character(),
                   pert = character(),
                   cell_line = character(), 
                   type = NULL,
                   where_clause = character(),
                   gold = FALSE,
                   inferred = TRUE,
                   fields = character(),
                   controls = FALSE,
                   cl = NULL,
                   verbose = FALSE) 
          {
            if (length(ids) && length(where_clause)) {
              stop("loadL1K function cannot be called with both ids and ",
                   "where_clause.")
            }
            if (!length(ids) && 
               !length(pert) && 
               !length(type) &&
               !length(cell_line) &&
               !length(where_clause) && 
               !gold) {
              if (verbose) message("Loading expression data...")
              if (inferred) {
                data <- readGCTX(x)
              } else {
                data <- readGCTX(x[seq_len(978),])
              }
            } else {
              if (length(where_clause)) {
                if (verbose)
                  message("Querying and loading expression data...")
                if (length(pert)) where_clause$pert_iname <- pert
                if (length(cell_line)) where_clause$cell_id <- cell_line
                if (length(type)) where_clause$pert_type <- type
                if (gold) where_clause$is_gold = TRUE
                ids <- clueInstances(x, where_clause = where_clause)
              }
              
              coln <- colnames(x)
              ix <- which(coln %in% ids)
              
              if (!length(ix)) {
                warning("Specified where_clause or ids returned no results. ",
                        "Please verify your query/ids.")
                return(NULL)
              }
              
              if (inferred) {
                data <- readGCTX(x[, ix])
              } else {
                data <- readGCTX(x[seq_len(978), ix])
              }
              if (ncol(data) != length(ids)) {
                message(length(ids), " instances requested, but only ",
                        ncol(data), " instances present in gctx file.")
              }
              rm(coln)
            }
            
            if (controls) {
              if (verbose) message("Loading control data...")
              ids <- controls(x, ids = base::colnames(data), 
                              verbose = verbose,
                              cl = cl)
              ix = which(colnames(x) %in% ids$distil_id)
              if (inferred) {
                ctrl <- readGCTX(x[, ix])
              } else {
                ctrl <- readGCTX(x[seq_len(978), ix])
              }
              data <- cbind(data, ctrl)
            }
            ix <- match(base::colnames(data), metadata(x)$inst_id)
            info <- metadata(x[,ix])
            rownames(info) <- info$inst_id
            info$distil_id <- info$inst_id
            
            if (length(fields)) {
              info <- info[, which(base::colnames(info) %in% fields)]
            }

            if (!all.equal(base::rownames(info), 
                           as.vector(base::colnames(data)))) {
              stop("Rownames of metadata and colnames of expression data ",
                   "did not match. Aborting.")
            }
            SummarizedExperiment::SummarizedExperiment(
              assays = list(exprs = data), 
              colData = info)
            
          })
