#' @export
Slinky$methods(colnames = function(file = NULL, index = NULL) {
    "Retrieve column names from LINCS gctx datafile
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{file} Path to gctx file.  May be omitted if specified when
            Slinky object is created.}
        \\item{\\code{index} Optional list providing which colnames to return.}
    }}
    \\subsection{Return Value}{Names of columns from gctx file}
    \\subsection{Details}{The gctx file is an HDF5 formatted file with several
        sections (groups) containing the column and row level metadata as well
        as the expression data itself.  If `index` is provided, it should be a
        list of extent one providing a vector of column indices for which to
        return ids.  For example, index=list(c(1,2,3,10)).}"

    .self$close()  # cleanup in case there was a bad exit previously
    if ((length(index) & !typeof(index) == "list") || length(index) > 1) {
        stop("Index must be a list of extent 1 provided indices of instanceIds",
            "to return")
    }

    if (!length(file) && !length(.self$.gctx)) {
        stop("You must specify path to gctx file.")
    }

    if (!length(file)) file = .self$.gctx

    if (!length(index)) {
        info <- rhdf5::h5dump(file, load = FALSE)
        index <- list(1:info$`0`$META$COL$id$dim)
    }

    rhdf5::h5read(file, name = "0/META/COL/id", index = index)
})

#' @export
Slinky$methods(rownames = function(file = NULL, index = NULL) {
    "Retrieve rown names from LINCS gctx datafile
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{file} Path to gctx file.  May be omitted if specified when
            Slinky object is created.}
        \\item{\\code{index} Optional list providing which rownames to return.}
    }}
    \\subsection{Return Value}{Names of rows from gctx file}
    \\subsection{Details}{The gctx file is an HDF5 formatted file with several
        sections (groups) containing the column and row level metadata as well
        as the expression data itself.  If `index` is provided, it should be a
        list of extent one providing a vector of row indices for which to
        return ids.  For example, index=list(c(1,2,3,10)).}"

    .self$close()  # cleanup in case there was a bad exit previously
    if ((length(index) & !typeof(index) == "list") || length(index) > 1) {
        stop("Index must be a list of extent 1 provided indices of rows to ",
            "return")
    }

    if (!length(file) && !length(.self$.gctx)) {
        stop("You must specify path to gctx file.")
    }

    if (!length(file)) file = .self$.gctx

    if (!length(index)) {
        info <- rhdf5::h5dump(file, load = FALSE)
        index <- list(1:info$`0`$META$ROW$id$dim)
    }
    rhdf5::h5read(file, name = "0/META/ROW/id", index = index)
})

#' @export
Slinky$methods(readGCTX = function(file = NULL, index = NULL) {
    "Read portions of data matrix from LINCS gctx datafile
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{file} Path to gctx file.  May be omitted if specified
            when Slinky object is created.}
        \\item{\\code{index} list of extent 2 providing which columns and rows
            to return.  Since the LINCS datasets are relatively massive, the
            index argument is required to avoid inadvertently trying to slurp
            the entire dataset into memory.}
    }}
    \\subsection{Return Value}{Matrix of expression data with rownames and
        colnames appropriately set}
    \\subsection{Details}{Failure to close connections will result in warning
        messages at garbage collection. Calling close prior to removing object
        will close connections and prevent these messages.}"

    .self$close()  # cleanup in case there was a bad exit previously

    if (!length(file) && !length(.self$.gctx)) {
        stop("You must specify path to gctx file.")
    }

    if (!length(index)) {
        stop("You must provide index argument.  E.g. index=list(1:10, 1:5)")
    }

    if ((length(index) & !typeof(index) == "list") || length(index) != 2) {
        stop("Index must be a list of extent 2 provided indices of ",
            "instanceIds to return, e.g. 'list(1:10, 1:50)'")
    }

    if (!length(file)) file = .self$.gctx

    data <- rhdf5::h5read(file, name = "0/DATA/0/matrix", index = index)
    colnames(data) <- as.vector(.self$colnames(file, index = list(index[[2]])))
    rownames(data) <- as.vector(.self$rownames(file, index = list(index[[1]])))
    data
})

#' @export
Slinky$methods(close = function() {
    "Close any open HDF5 (gctx) file connections.
    \\subsection{Return Value}{None.  Called for side effect of closing
        connections.}"
    rhdf5::H5close()
})
