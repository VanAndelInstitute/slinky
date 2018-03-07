#' @export
Slinky$methods(
    ks = function(data) {
        "Calculate KS based enrichment score
        \\subsection{Parameters}{
        \\itemize{
        \\item{\\code{ks} An N x M matrix of gene expression values (typically
                 robust zscores) for N genes and M samples. Rownames should
                 be the gene ids.}
        }}
        \\subsection{Return Value}{Vector of KS scores for each gene.}
        \\subsection{Details}{The the original CMAP paper by Lamb et al. used
            a KS random walk based statistic for calculating enrichment.  Here
            we use the same metric to identify differentially expressed genes.
            The question arises how to summarize accross replicates.
            We use a \"rank of ranks\" approach.  With this
            method, each sample is ranked, and then the entire matrix is
            ranked to create a single vector.  The KS statistic is then
            calculated for each gene based on the position of their replicate
            values in the vectorized matrix.  In this way, genes that
            consistently have high ranks within each sample will have a higher
            KS score than those that are inconsistently ranked.}"

    .ks <- function(ix, n) {
        scores <- -rep(1/(n-length(ix)), n)
        inc <- 1/length(ix)

        # need to account for ties
        ix <- floor(ix)
        scores[ix] <- 0
        for(i in ix) {
            scores[i] = scores[i] + inc
        }
        if(-min(cumsum(scores)) >= max(cumsum(scores))) {
            return(min(cumsum(scores)))
        } else {
            return(max(cumsum(scores)))
        }
    }

    genes <- base::rownames(data)

    # rank in descending order
    ranks <- apply(-data, 2, rank)
    ranks <- rank(as.vector(ranks))
    genes <- rep(base::rownames(data), base::ncol(data))

    tapply(ranks, INDEX = genes, FUN = .ks, length(ranks))
})
