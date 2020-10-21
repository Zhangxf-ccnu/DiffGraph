\name{Dtrace}
\alias{Dtrace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lasso penalized D-Trace loss
}
\description{
This function is implemented to identify gene network rewiring using the Dtrace model.
}
\usage{
Dtrace(X, lambda, covType = "pearson", tol = 1e-05, maxiter = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A list (length=2) of input matrices. They can be data matrices or covariance matrices. If every matrix in the X is a symmetric matrix, the matrices are assumed to be covariance matrices among genes. If the input are data matrices, the rows represent the observations, and the columns represents the genes.
}
  \item{lambda}{
A nonnegative number. This tuning paramter controls the sparsity level of the estimated differential network.
}
  \item{covType}{
A parameter to decide which approach is selected to compute the sample covariance matrices.

If covType = "pearson", it means that we compute (when input X represents data) the sample covariance matrices using Pearson correlation and that Dtrace is implemented based on Gaussian graphical models.  If data is normal, we suggest covType = "pearson".

If covType = "kendall" or "spearman", it means that we compute (when input X represents data) the sample covariance matrices using kendall's tau correlation or Spearman's rho correlation and that Dtrace is implemented based on nonparanormal graphical models. If data is non-normal, we suggest covType = "kendall" or "spearman".
}
  \item{tol}{
The tolerance parameter for convergence criteria.
}
  \item{maxiter}{
The maximum number of iterations.
}
}
\details{
Lasso penalized D-Trace loss (Dtrace) estimates the precision matrix difference directly without estimating the condition specific precision matrices. One major difference of Dtrace as compared to FGL is that the Dtrace method only requires sparsity assumption on the differential network while FGL applies sparsity penalties to both the condition specific networks and the differential network. Dtrace also has an advantage of using much smaller sample size to achieve competitive performance.

We solve the optimization problem of Dtrace using the accelerated proximal gradient method. 
}
\value{
\item{Delta}{The estimated difference of condition specific precision matrices.}
\item{Delta.graph.full}{The estimated differential network over all nodes.}
\item{Delta.graph.connected}{The estimated differential network over only the connected nodes.}
}

\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Xiaohua Hu and Hong Yan (2017), DiffGraph: An R package for identifying gene
network rewiring using differential graphical models.\cr

Huili Yuan, Ruibin Xi, Chong Chen and Minghua Deng (2017). Differential network analysis via lasso penalized D-trace loss. Biometrika, 104(4), 755–770. \cr

Dechao Tian, Quanquan Gu and Jian Ma (2016). Identifying gene regulatory network rewiring using latent
differential graphical models. Nucleic Acids Research, 44(17), e140.
}
\author{
\packageAuthor{DiffGraph}

Maintainer: \packageMaintainer{DiffGraph}
}



\seealso{
\code{\link{FGL}}, \code{\link{PNJGL}}, \code{\link{pDNA}}, \code{\link{TCGA.BRCA}}, \code{\link{TCGA.GBM}}
}


\examples{

# Identify differential network between breast cancer subtypes
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
dtrace.results= Dtrace(X, 0.45, covType = "spearman")
net.dtrace = dtrace.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.dtrace, vertex.size= degree(net.dtrace)*1.5, layout =layout_with_fr, 
    vertex.color="red", vertex.label.cex=0.8, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.dtrace <- tkplot.getcoords(tkid)
plot(net.dtrace, layout=l.dtrace, vertex.size= degree(net.dtrace)*1.5, vertex.color="red", 
    vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")



# Identify differential network between glioblastoma subtypes
data(TCGA.GBM)
X = TCGA.GBM$X[1,]
dtrace.results= Dtrace(X, 0.35, covType = "spearman")
net.dtrace = dtrace.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.dtrace, vertex.size= degree(net.dtrace)*1.5, layout =layout_with_fr,
    vertex.color="red", vertex.label.cex=0.8, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.dtrace <- tkplot.getcoords(tkid)
plot(net.dtrace, layout=l.dtrace,  vertex.size= degree(net.dtrace)*1.5,  vertex.color="red", 
    vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
}
