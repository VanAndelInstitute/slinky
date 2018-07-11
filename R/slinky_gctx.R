#' readGCTX
#' 
#' Read portions of data matrix from LINCS gctx datafile
#' 
#' @param x a Slinky Object
#' @return Matrix of expression data with rownames and
#'    colnames appropriately set. If a subset of the data is desired, 
#'    subset the slinky object itself, not the resulting data matrix. 
#'    That is, \code{data <- readGCTX(x[1:50,1:500])} will be 
#'    MUCH faster than \code{ndata <- readGCTX(x)[1:50, 1:500]}.
#' @name readGCTX
#' @rdname readGCTX
setGeneric("readGCTX",
           function(x) {
             standardGeneric("readGCTX")
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
#' data <- readGCTX(sl[1:20,1:5])
#' 
#' @rdname readGCTX
#' @exportMethod readGCTX
#' @importFrom rhdf5 h5closeAll
#' @importFrom rhdf5 h5read
#' @aliases readGCTX
setMethod("readGCTX", signature(x = "Slinky"),
          function(x) 
          {
            
            closeAll(x)  # cleanup in case there was a bad exit previously
            
            if (!length(x@gctx)) {
              stop("You must specify path to gctx file when ",
                    "creating Slinky object.")
            }
            
            data <- rhdf5::h5read(x@gctx, name = "0/DATA/0/matrix", 
                                  index = x@.index)
            closeAll(x)
            colnames(data) <- as.vector(colnames(x))
            rownames(data) <- as.vector(rownames(x))
            closeAll(x) # release lock
            data
          })

#' closeAll
#' 
#' Close any open HDF5 (gctx) file connections.
#' @param x a Slinky Object
#' @return None.  Called for side effect of closing
#'    connections.
#' @name closeAll
#' @rdname closeAll
setGeneric("closeAll", 
           function(x) { standardGeneric("closeAll")})
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
#' closeAll(sl)
#' @rdname closeAll
#' @exportMethod closeAll
#' @aliases close
setMethod("closeAll", signature = signature("Slinky"),
          function(x)
          {
            rhdf5::h5closeAll()
          }
)
