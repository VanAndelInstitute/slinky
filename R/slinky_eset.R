#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @export
Slinky$methods(toEset = function(gctx=NULL, index=NULL, where_clause=NULL, fields=NULL, info_file=NULL, verbose=FALSE) {
  "Convert data from gctx file to eset, pulling metadata from various sources
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{gctx} Path to gctx file.  May be omitted if already set when Slinky object was instantiated.}
  \\item{\\code{index} Index of extent 2 giving the row and column indices to pull from the gctx file.}
  \\item{\\code{where_clause} Query to use to determine which columns to pull from gctx file.}
  \\item{\\code{fields} Fields to include in the expression set's phenodata.  Default is all available.}
  \\item{\\code{info_file} Optional path to GSE92742_Broad_LINCS_inst_info.txt or other identically formatted metadata file.  
                           that covers the requested instances.  Default is to look in the library's 
                           installation folder and current working directory, and then download it if not found.}
  \\item{\\code{verbose} Do you want to know how things are going?  Default is FALSE}}}
  \\subsection{Return Value}{Object of type eSet containing expression and pheno data.}
  \\subsection{Details}{You must specify either index or where_clause to avoid inadvertently slurping in entire gctx file.}"
  
  if(!length(info_file)) {
    info_file = .self$.locateInfo()
  }
  
  if(!length(info_file)) stop ("Could not locate or download LINCS info file.  Consider retrying with info_file argument set.")

  if(verbose) cat("\n\nLoading metadata...")
  info <- read.delim(info_file, as.is=TRUE, header=TRUE)
  
  if(length(index)) {
    if(length(gctx)) {
      if(verbose) cat("\n\nLoading expression data...")
      data <- .self$readGCTX(gctx, index=index)
    } else {
      if(verbose) cat("\n\nLoading expression data...")
      data <- .self$readGCTX(index=index)
    }
  } else {
    if(!length(where_clause)) {
      stop("Either index or where clause must be specified.")
    }
    if(verbose) cat("\n\nQuerying and loading expression data...")
    inst <- .self$query.instances(where_clause=where_clause)
    coln <- .self$colnames(gctx)
    ix <- which(coln %in% inst)
    data <- .self$readGCTX(index=list(NULL, ix))
    rm(coln)
  }
  ix <- match(base::colnames(data), info$inst_id)
  info <- info[ix,]
  rownames(info) <- info$inst_id
  info$distil_id <- info$inst_id
  if(length(fields)) {
    info <- info[,which(base::colnames(info) %in% fields)]
  }
  if(!all.equal(base::rownames(info), as.vector(base::colnames(data)))) {
    stop("Rownames of metadata and colnames of expression data did not match.  Aborting.")
  }
  
  pd <- new("AnnotatedDataFrame", data=info)
  Biobase::ExpressionSet(assayData=data,
                phenoData=pd)
})

#' @export
Slinky$methods(.locateInfo = function(path=NULL, verbose=TRUE) {
  "Internal function to locate LINCS instance info datafile and download if missing.  
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{path} Optional path of where to look.  Default is to look in installation directory and current working directory.}
  \\item{\\code{verbose} Do you want to know how things are going?  Default is FALSE}}}
  \\subsection{Return Value}{Path to info file (downloaded if missing as a side effect)}
  \\subsection{Details}{Use the `download` function for general use.}"
  
  # this is the instance level metadata file from LINCS
  fn <- "GSE92742_Broad_LINCS_inst_info.txt"
  loc <- NULL
  
  # try to locate metadata file
  info_lib <- file.exists(system.file("extdata", fn, package="slinky"))
  info_local <- file.exists(paste(getwd(), fn, join="/"))
  if(!length(path)) {
    info_other = FALSE;
  } else {
    info_other <- (file.exists(path) && !file.info(path)$isdir) || file.exists(paste(path, fn, sep="/"))
  }
  
  # if cannot find it, download it.
  if(!info_lib && !info_local && !info_other) {
    if(verbose) message("\n\nInstance info file not found.  Downloading...\n\n")
    download("info", prompt=FALSE)
    
    # try to move it into pacakge installation directory for future reference
    dest <- paste(system.file("extdata", package="slinky"), fn, sep="/")
    tryCatch({
      if(verbose) cat("\nAttempting to move into lib installation directory.")
      file.rename(fn, dest)
      loc <- dest
    }, 
    error = function(e) {
      cat("Could not move file into lib directory.  Will use from current working directory instead.")
      loc <- fn
    })
  } else if (info_local) {
    loc <- paste(getwd(), fn, join="/")
  } else if (info_other) {
    if(file.exists(path) && !file.info(path)$isdir) {
      loc <- path;
    } else if (file.exists(paste(path, fn, sep="/"))) {
      loc <- paste(path, fn, sep="/")
    }
  } else if(info_lib) {
    loc <- system.file("extdata", fn, package="slinky")
  }
  return(loc[1])
})


# sl$toEset(where_clause=list("pert_desc"="LMNA", "pert_type"="trt_sh", "is_gold"=TRUE))
