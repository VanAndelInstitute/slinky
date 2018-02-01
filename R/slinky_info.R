#' @export
Slinky$methods(loadInfo = function(file=NULL, verbose=FALSE, reload=FALSE) {
  "Load the info file into memory if not already loaded.  
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{file} Optional complete path to file to load.  If omitted, default file from GEO will be located, or downloaded if absent.}
  \\item{\\code{verbose} Do you want to know how things are going?  Default is FALSE.}
  \\item{\\code{reload} If metadata is already loaded, should we reload it? Default is FALSE.}}}
  \\subsection{Return Value}{None.  Called for side effect of loading .self$metadata attribute.}
  \\subsection{Details}{After complete, metadata from info file will be in .self$metadata.}"
  if(length(.self$metadata) && !reload) return;
  
  if(!length(file)) {
    file = .self$.locateInfo(verbose = verbose);
  }
  
  tryCatch({
    if(verbose) message("\nAttempting to load info file.  This may take a couple minutes (but only needs to be done once per session).")
    .self$metadata = read.delim(file, as.is = TRUE);
    .self$.info = file;
    
  }, 
  error = function(e) {
    print(e)
    if(length(file) || .self$.info != "GSE92742_Broad_LINCS_inst_info.txt") {
      message(paste0("Could not load specified info file: ", file, ". Please verify file path."))
    } else {
      message(paste0("Could not load default info file.  Please verify internet connectivity."))
    }
  })
})


#' @export
Slinky$methods(.locateInfo = function(verbose=FALSE) {
  "Internal function to locate LINCS instance info datafile and download if missing.  
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{path} Optional path of where to look.  Default is to look in installation directory and current working directory.}
  \\item{\\code{verbose} Do you want to know how things are going?  Default is FALSE}}}
  \\subsection{Return Value}{Path to info file (downloaded if missing as a side effect)}
  \\subsection{Details}{Use the `download` function for general use.}"
  
  # this is the instance level metadata file from LINCS
  fn <- .self$.info # defaults to "GSE92742_Broad_LINCS_inst_info.txt"
  loc <- NULL
  
  # try to locate metadata file
  info_lib <- file.exists(system.file("extdata", fn, package="slinky"))
  info_wd <- file.exists(paste(getwd(), fn, join="/"))
  info_user <- file.exists(fn)
  # if cannot find it, download it.
  if(!info_lib && !info_wd && !info_user) {
    if(.self$.info != "GSE92742_Broad_LINCS_inst_info.txt") {
      stop("Could not locate specified info file.  \nSet info attribute of Slinky object to NULL to fetch default file from GEO.")
    }
    
    if(verbose) message("Instance info file not found.  Downloading from GEO...")
    download("info", prompt=FALSE)
    
    # try to move it into pacakge installation directory for future reference
    dest <- paste(system.file("extdata", package="slinky"), fn, sep="/")
    tryCatch({
      if(verbose) cat("\nAttempting to move into lib installation directory.")
      file.copy(fn, dest) # file.rename dies accross devices
      unlink(fn)
      loc <- dest
    }, 
    error = function(e) {
      cat("Could not move file into lib directory.  Will use from current working directory instead.")
      loc <- fn
    })
  } else if (info_wd) {
    loc <- paste(getwd(), fn, join="/")
  } else if (info_user) {
    loc <- fn;
  } else if(info_lib) {
    loc <- system.file("extdata", fn, package="slinky")
  }
  return(loc[1])
})

#' @export
#' @importFrom foreach %dopar% foreach
#' @importFrom dplyr group_by summarize
# TODO: This function is pretty slow due to repeated calls to clue.io.  Need to optimize.
Slinky$methods(controls = function(ids, verbose=FALSE, cl=NULL) {
  "Fetch the same plate control samples applicable for given ids (distil_id).  Expects that the specified ids have pert_type of trt_sh or trt_cp.
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{ids} The distil_id(s) to lookup.}
  \\item{\\code{cl} Optional cluster object to parallelize this operation. If verbose is TRUE, use this pattern in order for progress bar to update:
  \\code{cl <- parallel::makeCluster(4, outfile=\"\")}}
  \\item{\\code{verbose} Do you want to know how things are going? Default is FALSE.}
  }}
  \\subsection{Return Value}{The name of the vehicle control for the queried perturbagen(s).}
  \\subsection{Details}{For a given set of distil_ids, this function finds the distil_ids for the corresponding control samples based on the 
                        the pert_type and (for trt_cp) the specified vehicle.  The returned dataframe can be used, among other things, to create 
                        a control dataset for differential expression or other analysis.  See also de.by.plate.}"
  if(!length(.self$metadata)) {
    .self$loadInfo(verbose=verbose);
  }
  
  profs <- .self$clue("profiles", 
                      fields=c("pert_id", "pert_iname", "pert_vehicle", "pert_type", "rna_plate"), 
                      verbose=verbose, 
                      cl=cl, 
                      ids=ids)
  cond <- profs %>% group_by(rna_plate, pert_type, pert_vehicle) %>% summarize(unique(pert_type))
  if(verbose) {
    message("\nIdentifying control samples.")
    flush.console()
    pb <- txtProgressBar(min = 0, 
                         max = nrow(cond), 
                         initial = 0, 
                         width = 100, 
                         style = 3)
    i <- 0
  }
  if(length(cl)) {
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  ctrls <- foreach(i=1:nrow(cond), .combine=dplyr::bind_rows, .export = c("sl")) %dopar% {
    if(cond[i,]$pert_type == "trt_sh" || cond[i,]$pert_type == "trt_oe") {
      ix <- which(.self$metadata$pert_type == "ctl_vector" &
                  .self$metadata$rna_plate == cond[i,]$rna_plate)
    } else if(cond[i,]$pert_type == "trt_cp") {
      ix <- which(.self$metadata$pert_iname == cond[i,]$pert_vehicle &
                  .self$metadata$rna_plate == cond[i,]$rna_plate)
    } else {
      stop(paste0("Controls can not be retrieved for pert_type ", cond[i,]$pert_type, " (only trt_sh, trt_oe, or trt_cp supported)."))
    }
    if(verbose) {
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    prof = .self$metadata[ix,c(1,5,6,2)]
    base::colnames(prof)[1] <- "distil_id"
    return(prof)
  }
  ctrls
})


