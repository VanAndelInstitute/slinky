#' A Reference Class to encapsulate LINCS methods
#'
#' @field .user_key private field populated at instantiation
#' @field .base private field populated at instantiation
#' @field .gctx private field optionally populated at instantiation specifying
#'         path to gctx file
#' @import methods
#' @examples
#'
#' # for build/demo only.  You MUST use your own key when using the slinky
#' # package.
#' user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
#'                           as='parsed')$user_key
#' sl <- Slinky$new(user_key,
#'                  system.file('extdata', 'demo.gctx',
#'                       package='slinky'),
#'                  system.file('extdata', 'demo_inst_info.txt',
#'                      package='slinky'))
#' amox_gold <- sl$clue.instances(where_clause=list('pert_type'='trt_cp',
#'                  'pert_iname'='amoxicillin',
#'                  'cell_id' = 'MCF7',
#'                  'is_gold'=TRUE), poscon = 'omit')
#' amox_gold_sumexp <- sl$toSummarizedExperiment(ids = amox_gold)
#'
#' @exportClass Slinky
#' @export Slinky
Slinky <- setRefClass(
    "Slinky",
    fields = list(
        .user_key = "character",
        .base = "character",
        .gctx = "ANY",
        .info = "ANY",
        metadata = "ANY"
    ),
    methods = list()
)

#' @export
Slinky$methods(
    initialize = function(key = NULL,
                        gctx = NULL,
                        info = NULL) {
    "Create a Slinky object
    \\subsection{Parameters}{
    \\itemize{
    \\item{\\code{key} Your clue.io user key.  If not specified, CLUE_API_KEY
    environment variable must be set.}
    \\item{\\code{gctx} Optional path to gctx file containing data you want to
    work with.  Can be specified later if desired.}
    \\item{\\code{info} Optional path to info file containing metadata
    describing the instances in the gctz file.}
    }}
    \\subsection{Return Value}{Slinky object}"

        if (!length(key)) {
            key = Sys.getenv("CLUE_API_KEY")
            if (nchar(key) < 10) {
                stop(
                    "No user key provided.  Either provide user_key argument ",
                    "or set CLUE_API_KEY env variable"
                )
            }
        }
        .self$.user_key = key
        .self$.base = "https://api.clue.io"
        if (length(gctx)) {
            .self$.gctx = gctx
        } else {
            .self$.gctx = NULL
        }
        if (length(info)) {
            .self$.info = info
        } else {
            .self$.info = "GSE92742_Broad_LINCS_inst_info.txt.gz"
        }
        .self$metadata = NULL
    }
)


#' @export
Slinky$methods(
    distilToPlate = function(x) {
        "Extract the plate id from distil_id
        \\subsection{Parameters}{
        \\itemize{
            \\item{\\code{x} A vector of distil_ids, a dataframe with
                distil_id element, or an SummarizedExperiment whose colData
                contains distil_id.}
        }}
        \\subsection{Return Value}{Slinky object}"

        if (class(x) == "SummarizedExperiment") {
            if (!length(x$distil_id)) {
                stop(
                    "SummarizedExperiment passed to distilToPlate function, ",
                    "colData did not contain distil_id element."
                )
            }
            plates <- gsub("(_X\\d{1,1})_.*", "\\1", x$distil_id)

        } else if (class(x) == "data.frame") {
            if (length(x$distil_id)) {
                plates <- gsub("(_X\\d{1,1})_.*", "\\1", x$distil_id)
            } else if (length(x$inst_id)) {
                plates <- gsub("(_X\\d{1,1})_.*", "\\1", x$inst_id)
            } else {
                stop(
                    "Dataframe passed to distilToPlate function, but it",
                    "did not contain distil_id or inst_id element."
                )
            }
        } else {
            plates <- gsub("(_X\\d{1,1})_.*", "\\1", x)
        }
        return(plates)

    }
)


#' @export
Slinky$methods(
  load = function(pert, 
                  cell_line = NULL,
                  type = c("trt_cp", "trt_sh", "trt_oe"), 
                  gold = TRUE, 
                  controls = TRUE,
                  verbose = FALSE) {
    "High level function to load samples based on specified perturbation.
    \\subsection{Parameters}{
    \\itemize{
      \\item{\\code{pert} Name of the desired perturbagen.}
      \\item{\\code{pert} Optional vector of cell lines to restrict results 
                          to.}
      \\item{\\code{type} Type of perturbation (one of trt_cp, trt_sh, trt_oe).}
      \\item{\\code{gold} Should we restrict to 'gold' instances? Default is 
                          yes.}
      \\item{\\code{controls} Include same-plate controls?  Default is yes.}
      \\item{\\code{verbose} Do you want to know how things are going?
              Default is FALSE}}}
    \\subsection{Return Value}{Summarized Experiment}
    \\subsection{Details}{This is a convenience wrapper to the
        \\code{toSummarizedExperiment} method which provides more granular 
        control.  Note that for \\code{type}, the options are \\code{trt_cp}
        (the default) for compound (i.e. drug) treated samples, 
        \\code{trt_sh} for short hairpin treated samples, and 
        \\code{trt_oe} for samples treated with over expression constructs.}"
    
    type = match.arg(type)
    where_clause = list(pert_iname = pert, 
                        pert_type = type)
    if (gold)
      where_clause$is_gold = TRUE
    if (length(cell_line))
      if(length(cell_line) > 1) {
        where_clause$cell_id = list(inq = c(cell_line))
      } else {
        where_clause$cell_id = cell_line
      }
    
    sl$toSummarizedExperiment(where_clause = where_clause,
                              controls = controls)
  }
)
