Slinky$methods(score = function(data, method, ...) {
  "Provide a uniform interface to scoring functions.
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{data} Matrix of expression data to score. }
  \\item{\\code{method} Scoring method to use.  Only \\code{ks} is presently supported.}
  \\item{\\code{...} Additional arguments for \\code{method}.}
  }}
  \\subsection{Return Value}{Vectors of scores, one per column of data}"
  if(method == "ks") {
    return(.self$.ks(data, ...))
  } else if(method == "xsum") {
    return(.self$.xsum(data, ...))
  } else {
    return(NULL)
  }
})

Slinky$methods(.ks = function(data, up, down, plot=FALSE) {
  "Calculate KS based enrichment statistic (CMAP method, Lamb, et al.)
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{data} Matrix of expression data to score. }
  \\item{\\code{up} List of upregulated gene ids (must match rownames(data)).}
  \\item{\\code{down} List of downregulated gene ids (must match rownames(data)).}
  \\item{\\code{plot} Whether to plot GSEA style plot for each column (not implemented).}
  }}
  \\subsection{Return Value}{Vectors of scores, all in (0,1), one per column of data}"
  f <- function(ix, n) {
    sc <- rep(-1/(n-length(ix)), n)
    # account for ties
    ix <- round(ix)
    sc[ix] <- 0
    for(i in ix) {
      sc[i] <- sc[i] + 1/length(ix) 
    }
    sc <- cumsum(sc) 
    if(max(sc) > -min(sc)) {
      return(max(sc))
    } else {
      return(min(sc))
    }
  }
  data <- apply(-data, 2, rank)
  up <- which(base::rownames(data) %in% up)
  down <- which(base::rownames(data) %in% down)
  apply(data, 2, function(x) {
    u <- f(x[up], length(x))
    d <- f(x[down], length(x))
    if(sign(u) == sign(d)) {
      return(0)
    } else {
      return((0.5*(abs(u) + abs(d))) * sign(u))
    }
  })
})

Slinky$methods(.xsum = function(data, up, down, n=100) {
  "Calculate Extreme Sum stat (Agerwal et al.))
  \\subsection{Parameters}{
  \\itemize{
  \\item{\\code{data} Matrix of expression data to score. }
  \\item{\\code{up} List of upregulated gene ids (must match rownames(data)).}
  \\item{\\code{down} List of downregulated gene ids (must match rownames(data)).}
  \\item{\\code{n} How many genes to consider for scoring among up and down regulated genes in each column of the data matrix.}
  }}
  \\subsection{Return Value}{Vectors of scores, all in one per column of data}"
  
  up.ix <- which(base::rownames(data) %in% up)
  down.ix <- which(base::rownames(data) %in% down)
  f <- function(a) {
    a_r <- rank(a)
    changed <- a * (a_r > ( length(a_r) - n) | a_r < n)
    sum(changed[up.ix]) - sum(changed[down.ix], na.rm=TRUE)  
  }
  apply(data, 2, f)
})