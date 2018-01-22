#' @importFrom Biobase exprs
#' @export
Slinky$methods(chDir = function(treated, control) {
  "Convenience wrapper to calculate Characteristic Direction Unity Vector based on two datasets.
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{treated} Expression data for treated samples, as `data.frame`, `matrix`, or `ExpressionSet`}
  \\item{\\code{control} Expression data for control samples, as `data.frame`, `matrix`, or `ExpressionSet`}}}
  \\subsection{Return Value}{Column matrix of characteristic direction scores for each gene.}
  \\subsection{Details}{There are a few steps involved in getting data formatted for the `chdirAnalysis` function.  
                        This function takes care of that for you.}"
  
  if(class(treated) == "ExpressionSet") {
    treated <- exprs(treated)
  }

  if(class(control) == "ExpressionSet") {
    control <- exprs(control)
  }

  dat <- cbind(treated, control)
  gg <- base::rownames(dat)
  dat <- apply(dat, 2, as.numeric)
  dat <- as.data.frame(dat)
  cl <- rep("1", base::ncol(dat))
  cl[1:base::ncol(treated)] <- "2"
  
  dat <- cbind(gg, dat, stringsAsFactors=FALSE)
  
  GeoDE::chdirAnalysis(dat, factor(cl))$chdirprops$chdir[[1]]
})

# crossprod(rr1, rr2)/sqrt(crossprod(rr1) * crossprod(rr2))
