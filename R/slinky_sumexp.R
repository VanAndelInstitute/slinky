#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
Slinky$methods(toSummarizedExperiment = function(gctx = NULL,
                                index = NULL,
                                ids = NULL,
                                where_clause = NULL,
                                inferred = FALSE,
                                fields = NULL,
                                controls = FALSE,
                                info_file = NULL,
                                cl = NULL,
                                verbose = FALSE) {
    "Convert data from gctx file to SummarizedExperiment, pulling metadata 
    from various sources
    \\subsection{Parameters}{
    \\itemize{
        \\item{\\code{gctx} Path to gctx file.  May be omitted if already set
            when Slinky object was instantiated.}
        \\item{\\code{index} Index of extent 2 giving the row and column
            indices to pull from the gctx file. Exactly 1 of index, ids, or
            where_clause should be specified.}
        \\item{\\code{ids} distil_ids to include in the Expression Set.
            Exactly 1 of index, ids, or where_clause should be specified.}
        \\item{\\code{where_clause} Query to use to determine which columns to
            pull from gctx file. Exactly 1 of index, ids, or where_clause
            should be specified.}
        \\item{\\code{inferred} Should the inferred (non-landmark) genes be
            included in the analysis? Default is FALSE. Ignored if index is
            specified.}
        \\item{\\code{fields} Fields to include in the expression set's
            phenodata. Default is all available.}
        \\item{\\code{controls} Should same-plate controls be identified and
            included?  Default is FALSE.}
        \\item{\\code{info_file} Optional path to
            GSE92742_Broad_LINCS_inst_info.txt or other identically formatted
            metadata file that covers the requested instances.  Default is to
            look in the library's installation folder and current working
            directory, and then download it if not found.}
        \\item{\\code{verbose} Do you want to know how things are going?
            Default is FALSE}}}
    \\subsection{Return Value}{Object of type eSet containing expression and
            pheno data.}
    \\subsection{Details}{You must specify exactly one of index, ids or
            where_clause to avoid inadvertently slurping in entire gctx file.}"

    if (sum(!length(index), !length(ids), !length(where_clause)) != 2) {
        stop("toSummarizedExperiment function requires exactly one of index,",
            "ids, or where_clause.")
    }

    if (!length(info_file)) {
        info_file = .self$.locateInfo()
    }

    if (!length(info_file))
        stop("Could not locate or download LINCS info file. ",
            "Consider retrying with info_file argument set.")

    if (verbose) message("Loading metadata...")

    .self$loadInfo(info_file, verbose = verbose)

    if (length(index)) {
        if (class(index) != "list") {
            stop("The 'index' argument must be a list, e.g. list(1:10, 1:50).")
        }
        if (verbose) message("Loading expression data...")
        data <- .self$readGCTX(gctx, index = index)
        if (ncol(data) != length(index[[2]])) {
            message(length(index[[2]]), " instances requested, but only ",
                ncol(data), " instances present in gctx file.")
        }
    } else {
        if (length(where_clause)) {
            if (verbose)
                message("Querying and loading expression data...")
            ids <- .self$clue.instances(where_clause = where_clause)
        }

        coln <- .self$colnames(gctx)
        ix <- which(coln %in% ids)

        if (!length(ix)) {
            warning("Specified where_clause or ids returned no results. ",
                "Please verify your query/ids.")
            return(NULL)
        }

        if (inferred) {
            data <- .self$readGCTX(index = list(NULL, ix))
        } else {
            data <- .self$readGCTX(index = list(seq_len(978), ix))
        }
        if (ncol(data) != length(ids)) {
            message(length(ids), " instances requested, but only ",
                ncol(data), " instances present in gctx file.")
        }
        rm(coln)
    }

    if (controls) {
        if (verbose) message("Loading control data...")
        ids <- .self$controls(ids = base::colnames(data), verbose = verbose,
            cl = cl)
        if (inferred) {
            rows = NULL

        } else {
            rows = seq_len(978)
        }
        cols = which(.self$colnames(gctx) %in% ids$distil_id)
        ctrl <- .self$readGCTX(gctx, index = list(rows, cols))
        data <- cbind(data, ctrl)
    }
    ix <- match(base::colnames(data), .self$metadata$inst_id)
    info <- .self$metadata[ix, ]
    rownames(info) <- info$inst_id
    info$distil_id <- info$inst_id
    if (length(fields)) {
        info <- info[, which(base::colnames(info) %in% fields)]
    }
    if (!all.equal(base::rownames(info), as.vector(base::colnames(data)))) {
        stop("Rownames of metadata and colnames of expression data did not ",
            "match. Aborting.")
    }
    SummarizedExperiment::SummarizedExperiment(assays = list(exprs = data), 
                                               colData = info)
})
