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
#' pacakge.
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
#' amox_gold_eset <- sl$toEset(ids = amox_gold)
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
            .self$.info = "GSE92742_Broad_LINCS_inst_info.txt"
        }
        .self$metadata = NULL
    }
)


#' @export
Slinky$methods(
    distil.to.plate = function(x) {
        "Extract the plate id from distil_id
        \\subsection{Parameters}{
        \\itemize{
            \\item{\\code{x} A vector of distil_ids, a dataframe with
                distil_id element, or an eset whose pData(eset) contains
                distil_id.}
        }}
        \\subsection{Return Value}{Slinky object}"

        if (class(x) == "eset") {
            if (!length(pData(x)$distil_id)) {
                stop(
                    "ExpressionSet passed to distil.to.plate function, ",
                    "pheno data did not contain distil_id element."
                )
            }
            plates <- gsub("(_X\\d{1,1})_.*", "\\1", pData(x)$distil_id)

        } else if (class(x) == "data.frame") {
            if (length(x$distil_id)) {
                plates <- gsub("(_X\\d{1,1})_.*", "\\1", x$distil_id)
            } else if (length(x$inst_id)) {
                plates <- gsub("(_X\\d{1,1})_.*", "\\1", x$inst_id)
            } else {
                stop(
                    "Dataframe passed to distil.to.plate function, but it",
                    "did not contain distil_id or inst_id element."
                )
            }
        } else {
            plates <- gsub("(_X\\d{1,1})_.*", "\\1", x)
        }
        return(plates)

    }
)
