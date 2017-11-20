\name{DiffGraph-package}
\alias{DiffGraph-package}
\alias{DiffGraph}
\docType{package}
\title{
\packageTitle{DiffGraph}
}
\description{
DiffGraph is an R package for identifying gene network rewiring between two different conditions from gene expression data.
}
\details{
DiffGraph integrates four influential differential graphical models which identify gene network rewiring through estimating the difference of condition specific precision matrices: fused graphical lasso (FGL),  lasso penalized D-Trace loss (Dtrace), perturbed-node joint graphical lasso (PNJGL) and prior information-induced differential network analysis (pDNA). DiffGraph also provides three approaches to compute the condition specific covariance matrices (Pearson, Spearman and Kendall).  To fit normal data using Gaussian graphical models, we suggest using Pearson correlation. To fit non-normal data using nonparanormal graphical models, Spearman's rho or Kendall's tau correlation is suggested. DiffGraph outputs the estimated differential networks as igraph graphs, which can be visualized in both non-interactive and interactive manners using the igraph R package.
}
\author{
\packageAuthor{DiffGraph}

Maintainer: \packageMaintainer{DiffGraph}
}
\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Xiaohua Hu and Hong Yan (2017), DiffGraph: An R package for identifying gene
network rewiring using differential graphical models.}


\seealso{
 \code{\link{FGL}}, \code{\link{Dtrace}}, \code{\link{PNJGL}}, \code{\link{pDNA},} \code{\link{TCGA.BRCA}}, \code{\link{TCGA.GBM}}
}


