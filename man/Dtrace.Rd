\name{Dtrace}
\alias{Dtrace}
\title{
Lasso penalized D-Trace loss
}
\description{
This function is implemented to identify gene network rewiring using the Dtrace model.
}
\usage{
Dtrace(X, lambda, covType = "pearson", tol = 1e-05, maxiter = 500, 
    rho = 0.1, rho.incr = 1.05, rho.max = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
   A list of input matrices. They can be data matrices or covariance matrices. If every matrix in the X is a symmetric matrix, the matrices are assumed to be covariance matrices among genes. If the input are data matrices, the rows represent the observations, and
   the columns represents the genes.
}
  \item{lambda}{
A nonnegative number. The hyperparameter controls the sparsity level of the estimated differential network.
}
  \item{covType}{
A parameter to decide which approach we choose to compute the sample covariance matrices.

If covType = "pearson", it means that we compute (when input X represents data directly) the sample covariance matrices using Pearson correlation and that Dtrace is implemented based on Gaussian graphical models.  If data is normal, we suggest covType = "pearson".

If covType = "kendall" or "spearman", it means that we compute (when input X represents data directly) the sample covariance matrices using kendall's tau correlation or Spearman's rho correlation and that Dtrace is implemented based on nonparanormal graphical models. If data is non-normal, we suggest covType = "kendall" or "spearman".
}
  \item{tol}{
The tolerance parameter for convergence criteria.
}
  \item{maxiter}{
The maximum number of iterations for the ADMM algorithm.
}
  \item{rho}{
The penalty parameter in the ADMM algorithm.
}
  \item{rho.incr}{
The increase step parameter for varying penalty parameter rho.
}
  \item{rho.max}{
The maximum value of rho.
}
}
\details{Dtrace estimates the precision matrix difference directly without estimating the condition specific precision matrices. It contains a D-trace loss and a lasso penalty:
\deqn{
\min\limits_{\Delta = \Delta^T}   \  \left\{ \frac{1}{2} \left\langle   \Delta^2 ,  \frac{ S^{\left(1\right)} S^{\left(2\right)} + S^{\left(2\right)}S^{\left(1\right)} } {2 } \right\rangle - \left\langle \Delta, S^{\left(1\right)} -  S^{\left(2\right)} \right\rangle \right\} + \lambda \left\| \Delta \right\|_1,
}
where the first term is the D-trace loss function, the second term is the lasso penalty function, and \eqn{\lambda} is a nonnegative tuning parameter. Unlike FGL, Dtrace does not require assumption on the condition specific precision matrices. In addition, it has an advantage of using much smaller sample size to achieve competitive performance. The estimated difference of condition specific precision matrices \eqn{\hat{\Delta}} is used to construct the differential network.\cr

We solve the optimization problem using the ADMM algorithm. We accelerate the ADMM iterations by adaptively changing \eqn{\rho} in iterations. 
}
\value{
\item{Delta}{The estimated difference of condition specific precision matrices.}
\item{Delta.graph.full}{The estimated differential network over all nodes.}
\item{Delta.graph.connected}{The estimated differential network over only the connected nodes.}


}
\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Xiaohua Hu and Hong Yan (2017), DiffGraph: An R package for identifying gene
network rewiring using differential graphical models.\cr

Huili Yuan, Ruibin Xi, Chong Chen and Minghua Deng (2015). Differential network analysis via the lasso penalized d-trace loss. arXiv preprint arXiv:1511.09188.\cr

Dechao Tian, Quanquan Gu and Jian Ma (2016). Identifying gene regulatory network rewiring using latent
differential graphical models. Nucleic Acids Research, 44(17), e140
}
\author{
\packageAuthor{DiffGraph}

Maintainer: \packageMaintainer{DiffGraph}
}


\seealso{
\code{\link{FGL}}, \code{\link{PNJGL}}, \code{\link{pDNA}}, \code{\link{TCGA.BRCA}},\code{\link{TCGA.GBM}}
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
dtrace.results= Dtrace(X, 0.32, covType = "spearman")
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


