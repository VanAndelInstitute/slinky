#' An S4 class for working with LINCS data
#'
#' The \code{Slinky} class encapsulates details about location of the LINCS L1000 
#' data files as well as access credentials for the clue.io API (if desired).
#' It provides methods for querying and loading data from these resources.
#' The helper function \code{\link{Slinky}} is a simpler way to construct 
#' an object of this class
#' 
#' @slot  .index internal slot for mapping object to file data
#' @slot  base Base url for clue.io API.  
#' @slot  gctx gctx containing expression data (optional)
#' @slot  info info file containing metadata (optional)
#' @slot  metadata internal slot for storing metadata from info file, mapped 
#'        to gctx file and current index.
#' @slot  user_key clue.io API key (required unless CLUE_API_KEY env variable 
#'          is set)
#' @import methods
#' @import jsonlite
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
#' amox_gold <- clueInstances(sl, where_clause = list('pert_type' = 'trt_cp',
#'                  'pert_iname' = 'amoxicillin',
#'                  'cell_id' = 'MCF7',
#'                  'is_gold' = TRUE), poscon = 'omit')
#' amox_gold_sumexp <- loadL1K(sl, ids = amox_gold)
#'
#' @name Slinky-class
#' @exportClass Slinky
setClass("Slinky",
  representation(
    user_key = "character",
    gctx =  "character",
    info =  "character",
    metadata = "data.frame",
    base = "character",
    .index = "ANY"
  ),
  prototype(user_key = NA_character_, 
            gctx = NA_character_,
            info = NA_character_,
            metadata = data.frame(),
            base = NA_character_,
            .index = NULL),
  validity = function(object) {
  }
)
 
#' The \code{Slinky()} constructor returns a slinky object with defaults 
#' set where required. To access the clue.io API, you must either provide 
#' your key as the value for the arguemtn \code{user_key}, or set the 
#' CLUE_API_KEY environment variable prior to starting your R session.  
#' If no clue.io API access is required, simply specify \code{user_key = ""}. 
#' If \code{info} is not specified, it will be assumed that the file 
#' GSE92742_Broad_LINCS_inst_info.txt.gz is present in the working directory.
#' @param  user_key clue.io API key
#' @param  info info file containing metadata (optional)
#' @param  gctx gctx containing expression data (optional)
#' @return A Slinky object.
#' @rdname Slinky-class
#' @name Slinky
#' @export
#' @importFrom readr read_delim
#' @importFrom stats mad median
#' @importFrom utils flush.console read.delim txtProgressBar
Slinky = function(user_key = character(), 
                  gctx = character(), 
                  info = character()) {
  if (length(user_key) != 1) {
    key = Sys.getenv("CLUE_API_KEY")
    if (nchar(key) < 10) {
      stop(
        "No user key provided.  Either provide user_key argument ",
        "or set CLUE_API_KEY env variable"
      )
    }
    user_key <- key
  }
  
  if (!length(info)) {
    info <- "GSE92742_Broad_LINCS_inst_info.txt.gz"
  }
  
  meta <- rhdf5::h5dump(gctx, load = FALSE)
  ncol <- meta$`0`$META$COL$id$dim
  nrow <- meta$`0`$META$ROW$id$dim
  index <- list(seq_len(nrow), seq_len(ncol))
  
  base <- "https://api.clue.io"
  x <- new("Slinky", 
      user_key = user_key, 
      info = info, 
      base = base, 
      gctx = gctx,
      .index = index)       
  x@metadata <- .loadInfo(x)
  x
}



#' Slinky object dimensions
#' 
#' Get the number of rows and columns in L1000 data represented by 
#' Slinky object.
#' @return number of rows or columns of current (possibly subsetted) L1000 
#'     data set. 
#' @param x an object of class Slinky 
#' @name nrow
#' @rdname nrow
#' @exportMethod nrow
setGeneric("nrow")
#' @rdname nrow
#' @aliases nrow,ANY-method
setMethod("nrow", c("Slinky"), function(x) { return(length(x@.index[[1]]))})

#' @name ncol
#' @rdname nrow
#' @exportMethod ncol
setGeneric("ncol")
#' @rdname nrow
#' @aliases ncol,ANY-method
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
#' ncol(sl)
#' nrow(sl)
setMethod("ncol", c("Slinky"), function(x) { return(length(x@.index[[2]]))})

#' colnames
#' 
#' Retrieve column names from LINCS gctx datafile
#' @param x a Slinky Object
#' @param do.NULL Ignored (see \code{?base::colnames})
#' @param prefix Ignored (see \code{?base::colnames})
#' @return Names of columns from gctx file
#'          The gctx file is an HDF5 formatted file with several
#'           sections (groups) containing the column and row level metadata as well
#'           as the expression data itself.  Note that for best performance, if a 
#'           subset of colnames is desired, subset the slinky object itself, not
#'           the colnames, to avoid loading the entire set of colnames from the the 
#'           gctx file.  That is, \code{names <- colnames(x[,1:50])} will be 
#'           considerably faster than \code{names <- colnames(x)[1:50]}.
#'           The \code{do.NULL} and \code{prefix} arguments from 
#'           \code{base::colnames} do not apply here (as the slinky object will 
#'           always have column names), and will be silently ignored if provided.
#' @name colnames
#' @rdname colnames
setGeneric("colnames")
#' @exportMethod colnames
#' @rdname colnames
#' @aliases colnames,ANY-method
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
#' colnames(sl[,1:5])
setMethod("colnames", signature(x = "Slinky"),
          function(x) 
          {
            closeAll(x)  # cleanup in case there was a bad exit previously
            if (!length(x@gctx)) {
              stop("You must specify path to gctx file when creating ", 
                   "Slinky object.")
            }
            
            rhdf5::h5read(x@gctx, 
                          name = "0/META/COL/id", 
                          index = list(x@.index[[2]]))
          })

#' rownames
#' 
#' Retrieve row names from LINCS gctx datafile
#' @param x a Slinky Object
#' @param do.NULL Ignored (see \code{?base::rownames})
#' @param prefix Ignored (see \code{?base::rownames})
#' @return Names of rows from gctx file
#'           The gctx file is an HDF5 formatted file with several
#'           sections (groups) containing the column and row level metadata as well
#'           as the expression data itself.  Note that for best performance, if a 
#'           subset of rownames is desired, subset the slinky object itself, not
#'           the rownames, to avoid loading the entire set of rownames from the the 
#'           gctx file.  That is, \code{names <- rownames(x[,1:50])} will be 
#'           faster than \code{names <- rownames(x)[1:50]}.
#'           The \code{do.NULL} and \code{prefix} arguments from 
#'           \code{base::rownames} do not apply here (as the slinky object will 
#'           always have row names), and will be silently ignored if provided.
#' @name rownames
#' @rdname rownames
setGeneric("rownames")
#' @rdname rownames
#' @exportMethod rownames
#' @aliases rownames,Slinky-method
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
#' rownames(sl[1:5,])
setMethod("rownames", signature(x = "Slinky"),
          function(x) 
          {
            closeAll(x)  # cleanup in case there was a bad exit previously
            rhdf5::h5read(x@gctx, 
                          name = "0/META/ROW/id", 
                          index = list(x@.index[[1]]))
          })


#' Subsetting Slinky objects
#' 
#' @param x A Slinky Object
#' @param i row index
#' @param j column indes
#' Subsets a Slinky object.  This does not touch the data on file, it simply
#' adjusts the index slots in the resulting Slinky object to speed up 
#' subsequent data operations.
#' @rdname subset
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
#' colnames(sl[,1:5])
#' rownames(sl[1:5,1:5])
#' @return The subsetted Slinky object
setMethod("[", c("Slinky", "numeric", "numeric"),
          function(x, i, j) {
            x@.index = list(x@.index[[1]][i], 
                            x@.index[[2]][j])
            if (length(x@metadata))
              x@metadata <- x@metadata[j,]
            if (!.check_metadata(x)) {
              stop("Unexpected metadata mismatch.")
            }
            x
          })

#' @rdname subset
setMethod("[", c("Slinky", "missing", "numeric"),
          function(x, i, j) {
            x@.index = list(x@.index[[1]], 
                            x@.index[[2]][j])
            if (length(x@metadata))
              x@metadata <- x@metadata[j,]
            if (!.check_metadata(x)) {
              stop("Unexpected metadata mismatch.")
            }
            x
          })

#' @rdname subset
setMethod("[", c("Slinky", "missing", "missing"),
          function(x, i, j) {
            x
          })

#' @rdname subset
setMethod("[", c("Slinky", "numeric", "missing"),
          function(x, i, j) {
            x@.index = list(x@.index[[1]][i], 
                            x@.index[[2]])
            x
          })

