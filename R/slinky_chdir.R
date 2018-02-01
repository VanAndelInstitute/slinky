#' @importFrom Biobase exprs
#' @export
Slinky$methods(chDir = function(treated, control) {
  "Convenience wrapper to calculate Characteristic Direction Unity Vector based on two datasets.
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{treated} Expression data for treated samples, as `data.frame`, `matrix`, or `ExpressionSet`}
  \\item{\\code{control} Expression data for control samples, as `data.frame`, `matrix`, or `ExpressionSet`}
  \\item{\\code{plot} Do you want the projection and Top genes plots displayed? Default is FALSE.}
  }}
  \\subsection{Return Value}{Column matrix of characteristic direction scores for each gene.}
  \\subsection{Details}{There are a few steps involved in getting the data formatted for the `chdirAnalysis` function.  
                        This function takes care of that for you.  Also, the chdirSig function is not exported from 
                        the package GeoDE so we copy it here to be able to circumvent plotting (which may not be 
                        desired for the high throughput applications targeted by this package).}"
  
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
  #res <- GeoDE::chdirAnalysis(dat, factor(cl))$chdirprops$chdir[[1]]
  chdirSig(dat, factor(cl))$chdir[[1]]
})

# crossprod(rr1, rr2)/sqrt(crossprod(rr1) * crossprod(rr2))

# I want to avoid plotting for high throughput and/or headless analysis, but the workhorse 
# function called by GeoDE::chdirAnalysis is not exported.  So I reproduce it here.
#
# Copied verbatim from https://raw.githubusercontent.com/cran/GeoDE/master/R/ChDir-06.R
# Author: Neil R. Clark and Avi Ma'ayan
# Maintainer: Neil Clark <neil.clark@mssm.edu>
# Description: Given expression data this package calculate a multivariate
#    geometrical characterization of the differential expression and can also
#    perform gene-set enrichment.
# License: GPL-2

chdirSig <- function(data,sampleclass,gammas=list(1.0),nnull=10)
{
  
  pca1 <- prcomp(t(as.matrix(data[-1])))
  
  meanvec <- rowMeans(as.matrix(data[-1][sampleclass==2]))-rowMeans(as.matrix(data[-1][sampleclass==1]))
  
  n1 <- sum(sampleclass==1)
  n2 <- sum(sampleclass==2)
  
  cumsum <- pca1$sdev^2/sum(pca1$sdev^2)
  keepPC <- length(cumsum[cumsum>0.001])
  
  V <- pca1$rotation[,1:keepPC]
  R <- pca1$x[,1:keepPC]
  
  Dd <- (t(R[sampleclass==1,])%*%R[sampleclass==1,]+t(R[sampleclass==2,])%*%R[sampleclass==2,])/(n1+n2-2)
  
  sigma <- mean(diag(Dd))
  
  ShrunkMats <- lapply(gammas, function(x) solve(x*Dd + sigma*(1-x)*diag(keepPC)))
  
  b <- lapply(ShrunkMats, function(x) matrix(V%*%x%*%t(V)%*%meanvec,dimnames=list(c(as.list(as.character(data[[1]]))), 1)))
  
  b <- lapply(b, function(x) x/sqrt(sum(x^2)))
  
  b2dscale <- colMeans(R[sampleclass==2,1:2])- colMeans(R[sampleclass==1,1:2])
  
  b2dscale <- sqrt(sum(b2dscale^2))

  projchdir2d <-lapply(b, function(x) list( b2dscale*as.numeric(as.vector(x)%*%as.vector(V[,1])), b2dscale*as.numeric(as.vector(x)%*%as.vector(V[,2]))))  
  
  list(chdir=b,pca2d=R[,1:2],chdir_pca2d=projchdir2d) 
  
}
