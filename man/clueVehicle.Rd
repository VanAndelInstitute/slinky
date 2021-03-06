% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slinky_clue.R
\docType{methods}
\name{clueVehicle}
\alias{clueVehicle}
\alias{clueVehicle,Slinky-method}
\title{clueVehicle}
\usage{
clueVehicle(x, ids, verbose = FALSE)

\S4method{clueVehicle}{Slinky}(x, ids, verbose = FALSE)
}
\arguments{
\item{x}{a Slinky Object}

\item{ids}{The distil_id(s) to lookup.}

\item{verbose}{Do you want to know how things are going?
Default is FALSE.}
}
\value{
The name of the vehicle control for the queried
   perturbagen(s).
   This is a convenience wrapper to the profiles API
   which queries clue.io and unwraps response.
}
\description{
Fetch the vehicle control applicable to given ids (distil_id).  Expects
   that perturbagen is of type trt_cp.
}
\examples{
# for build/demo only.  You MUST use your own key when using the slinky
# package.
user_key <- httr::content(httr::GET('https://api.clue.io/temp_api_key'),
                          as='parsed')$user_key
sl <- Slinky(user_key,
                 system.file('extdata', 'demo.gctx',
                      package='slinky'),
                 system.file('extdata', 'demo_inst_info.txt',
                     package = 'slinky'))
amox <- clue(sl, where_clause = list("pert_iname" = "amoxicillin", 
                                    "cell_id" = "MCF7",
                                    "is_gold" = TRUE))
amox.ctrl <- clueVehicle(sl, amox$distil_id)

}
