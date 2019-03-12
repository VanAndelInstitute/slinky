#' metadata
#' 
#' The accessor function retrieves metadata from Slinky object. 
#' 
#' @param x a Slinky object
#' @return The accessor function returns a \code{data.frame} 
#' containing the metadata.
#'
#' @name metadata
#' @rdname metadata
setGeneric("metadata", function(x) standardGeneric("metadata"))
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
#' md <- metadata(sl[, 1:10])
#' 
#' @rdname metadata
#' @exportMethod metadata
#' @aliases metadata,Slinky-method
#' @rdname metadata
setMethod("metadata", signature(x = "Slinky"),
          function(x) 
          {
            return(x@metadata)
          })            


#' get_metadata
#' 
#' The accessor function retrieves metadata from Slinky object. 
#' 
#' As it turns out `metadata` was a poor choice for the accessor 
#' function because it can be masked if the user loads the 
#' `SummarizedExperiment` package after slinky.  So this provides 
#' an alternative.  Eventually `slinky::metadata` should be deprecated.
#' 
#' @param x a Slinky object
#' @return The accessor function returns a \code{data.frame} 
#' containing the metadata.
#'
#' @name get_metadata
#' @rdname get_metadata
setGeneric("get_metadata", function(x) standardGeneric("get_metadata"))
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
#' md <- get_metadata(sl[, 1:10])
#' 
#' @rdname get_metadata
#' @exportMethod get_metadata
#' @aliases get_metadata,Slinky-method
#' @rdname get_metadata
setMethod("get_metadata", signature(x = "Slinky"),
          function(x) 
          {
            return(x@metadata)
          })            



#' Non-exported functions from the slinky package.
#' 
#' .loadInfo is an internal function to load the info file into memory if 
#' not already loaded. 
#' Metadata from info file will then be available in \code{@metadata} slot. 
#' @param x  An object of class Slinky
#' @param file  Optional complete path to file to load.  If
#'    omitted, default file from GEO will be located, or downloaded
#'    if absent.
#' @param verbose  Do you want to know how things are going?
#'    Default is FALSE.
#' @param reload  If metadata is already loaded, should we reload
#'    it? Default is FALSE.
#' @return None.  Called for side effect of loading
#'    metadata slot.
#'
#' @name loadInfo
#' @aliases .loadInfo
#' @rdname slinky-internal
#' @noRd
setGeneric(".loadInfo",
           function(x, 
                    file = NULL, 
                    verbose = FALSE,
                    reload = FALSE) {
             standardGeneric(".loadInfo")
           }
)
#' @seealso metadata
#' @rdname slinky-internal
#' @aliases loadInfo, Slinky-method
#' @noRd
setMethod(".loadInfo", signature(x = "Slinky"),
          function(x, 
                   file = NULL, 
                   verbose = FALSE,
                   reload = FALSE) 
          {
            
            if (length(x@metadata) && !reload) return()
            
            if (!length(file)) {
              file = .locateInfo(x)
            }
            tryCatch({
              if (verbose)
                message("\nAttempting to load info file.  This may take a ",
                        "minute (but only needs to be done once per session).")
              md <- read.delim(file, as.is = TRUE)
            }, error = function(e) {
              print(e)
              if (length(file) ||
                  !grepl("GSE92742_Broad_LINCS_inst_info.txt.gz", x@info)) {
                message("Could not load specified info file: ", file,
                        ". Please verify file path.")
              } else {
                message("Could not load default info file.  Please verify ",
                        "internet connectivity.")
              }
            })
            cn <- colnames(x)
            ids <- md$inst_id
            if (grepl(":.+:", cn[1])) { # level5
              gctx.trt <- gsub(".*?:(.*?-.*?)-.*", "\\1", cn)
              gctx.plate <- gsub(":.*", "", cn)
              
              info.plate <- gsub("_X.*", "", md$rna_plate)
              info.trt <- md$pert_id
              
              ix <- match(paste(gctx.trt, gctx.plate), paste(info.trt, info.plate))
              md <- md[ix,]
              md$inst_id <- cn
            } else {
              ix <- match(cn, ids)
              md <- md[ix,]
              if (!all.equal(as.character(cn), as.character(md$inst_id)))
                stop("inst_ids in info file did not match colnames of gctx ",
                     "file.  \nPlease ensure info file containes metadata for ", 
                     "all columns of gctx file.")
            }
            md
          })


#' .locateInfo is an Internal function to locate LINCS instance info datafile 
#' and download if missing.
#' Use the \code{\link{download}} to download various L1000 associated 
#' files, or the \code{\link{metadata}} accessor function which will 
#' automatically download the info file necessary.
#' @name .locateInfo
#' @rdname slinky-internal
#' @noRd
setGeneric(".locateInfo",
           function(x, 
                    verbose = FALSE){
             standardGeneric(".locateInfo")
           }
)
#' @rdname slinky-internal
#' @aliases locateInfo
#' @noRd
setMethod(".locateInfo", signature(x = "Slinky"),
          function(x, verbose = FALSE) 
          {
            # this is the instance level metadata file from LINCS
            fn <- x@info  # defaults to 'GSE92742_Broad_LINCS_inst_info.txt.gz'
            loc <- NULL
            
            # try to locate metadata file
            info_lib <- file.exists(system.file("extdata", 
                                                fn, 
                                                package = "slinky"))
            info_wd <- file.exists(paste(getwd(), fn, join = "/"))
            info_user <- file.exists(fn)
            
            # if cannot find it, download it.
            if (!info_lib && !info_wd && !info_user) {
              if (x@info != "GSE92742_Broad_LINCS_inst_info.txt.gz") {
                stop("Could not locate specified info file.  \nSet info 
                     attribute of Slinky object to NULL to fetch default 
                     file from GEO.")
              }
              
              if (verbose)
                message("Instance info file not found. Downloading from GEO.")
              download(x, "info", prompt = FALSE)
              
              # try to move it into pacakge installation directory for
              # future reference
              dest <- paste(system.file("extdata", package = "slinky"),
                            fn, sep = "/")
              tryCatch({
                if (verbose)
                  cat("\nAttempting to move into lib installation directory.")
                file.copy(fn, dest)  # file.rename dies accross devices
                unlink(fn)
                loc <- dest
              }, error = function(e) {
                message("Could not move file into lib directory.",   
                        "Will use from current working directory instead.")
                loc <- fn
              })
            } else if (info_wd) {
              loc <- paste(getwd(), fn, join = "/")
            } else if (info_user) {
              loc <- fn
            } else if (info_lib) {
              loc <- system.file("extdata", fn, package = "slinky")
            }
            x@info <- loc[1]
            return(loc[1])
          })


#' .check_metadata is a security blanket to make sure metadata remains in 
#' sync with expression data represented by Slinky object.
#' @return Logical indicating whether we are in sync.
#'    #' @name .check_metadata
#' @rdname slinky-internal
#' @noRd
setGeneric(".check_metadata",
           function(x, verbose = FALSE) {
             standardGeneric(".check_metadata")
           }
)
#' @rdname slinky-internal
#' @aliases .check_metadata
#' @noRd
setMethod(".check_metadata", signature(x = "Slinky"),
          function(x) 
          {
            return(all.equal(as.character(colnames(x)), 
                              as.character(metadata(x)$inst_id)))
          })
            
            
#' controls
#'    #' Fetch the same plate control samples applicable for given ids (distil_id).
#'    Expects that the specified ids have pert_type of trt_sh or trt_cp.
#' @param x  A slinky object
#' @param ids  The distil_id(s) to lookup.
#' @param cl  Optional cluster object to parallelize this
#'    operation. If verbose is TRUE, use this pattern in order for
#'    progress bar to update:
#'    \code{cl <- parallel::makeCluster(4, outfile=\"\")}
#' @param verbose  Do you want to know how things are going?
#'    Default is FALSE.
#' @return The name of the vehicle control for the queried
#'    perturbagen(s).
#'    \\For a given set of distil_ids, this function finds
#'    the distil_ids for the corresponding control samples based on the the
#'    pert_type and (for trt_cp) the specified vehicle.  The returned
#'    dataframe can be used, among other things, to create a control dataset
#'    for differential expression or other analysis.  See also diffexp.
#' @name controls
#' @rdname controls
setGeneric("controls",
           function(x, ids, verbose = FALSE, cl = NULL) {
             standardGeneric("controls")
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
#' amox_gold <- clueInstances(sl, where_clause = list("pert_type" = "trt_cp",
#'    "pert_iname" = "amoxicillin",
#'    "cell_id" = "MCF7",
#'    "is_gold" = TRUE), 
#'    poscon = "omit")
#' colnames(sl[,1:5])
#' rownames(sl[1:5,1:5])
#' ids.ctrl <- controls(sl, ids = amox_gold)$distil_id
#' 
#' @rdname controls
#' @exportMethod controls
#' @aliases controls,Slinky-method
#' @import foreach
#' @importFrom dplyr group_by summarize %>%
setMethod("controls", signature(x = "Slinky"),
          function(x, ids, verbose = FALSE, cl = NULL) 
          {
            
            profs <- clue(x, "profiles", 
                          fields = c("pert_id", "pert_iname",
                                     "pert_vehicle", "pert_type", 
                                     "rna_plate"),
                          verbose = verbose,
                          cl = cl, 
                          ids = ids)
            
            pert_type <- NULL # prevent no visible binding on R CMD CHECK
            rna_plate <- NULL # prevent no visible binding on R CMD CHECK
            pert_vehicle <- NULL # prevent no visible binding on R CMD CHECK
            cond <- profs %>% 
              group_by(rna_plate, pert_type, pert_vehicle) %>%
              summarize(unique(pert_type))
            
            if (verbose) {
              message("\nIdentifying control samples.")
              flush.console()
              pb <- txtProgressBar(min = 0, max = nrow(cond), initial = 0,
                                   width = 100, style = 3)
              i <- 0
            }
            
            if (length(cl)) {
              doParallel::registerDoParallel(cl)
            } else {
              foreach::registerDoSEQ()
            }
            
            ctrls <- foreach(i = seq_len(nrow(cond)), .combine = dplyr::bind_rows,
                             .export = c("sl")) %dopar% 
              {
               if (cond[i, ]$pert_type == "trt_sh" ||
                   cond[i, ]$pert_type == "trt_oe") {
                 ix <- which(x@metadata$pert_type == "ctl_vector" &
                               x@metadata$rna_plate == cond[i, ]$rna_plate)
                 
               } else if (cond[i, ]$pert_type == "trt_cp") {
                 ix <- which(x@metadata$pert_iname == cond[i,]$pert_vehicle &
                               x@metadata$rna_plate == cond[i, ]$rna_plate)
                 
               } else {
                 stop("Controls can not be retrieved for pert_type ",
                      cond[i, ]$pert_type,
                      " (only trt_sh, trt_oe, or trt_cp supported).")
               }
               
               if (verbose) {
                 i <- i + 1
                 utils::setTxtProgressBar(pb, i)
               }
               
               prof = x@metadata[ix, c(1, 5, 6, 2)]
               base::colnames(prof)[1] <- "distil_id"
               return(prof)
             }
            ctrls
          }
  )
