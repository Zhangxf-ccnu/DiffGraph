\name{PNJGL}
\alias{PNJGL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perturbed-node joint graphical lasso
}
\description{
This function is implemented to identify gene network rewiring using the PNJGL model.
}
\usage{
PNJGL(X, lambda1, lambda2, covType = "pearson", weights = "equal", 
    penalize.diagonal = FALSE, tol = 1e-05, maxiter = 500, rho = 0.1, 
    rho.incr = 1.05, rho.max = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A list (length=2)  of input matrices. They can be data matrices or covariance matrices. If every matrix in the X is a symmetric matrix, the matrices are assumed to be covariance matrices among genes. If the input are data matrices, the rows represent the observations, and
   the columns represents the genes.
}
  \item{lambda1}{
A nonnegative number. The hyperparameter controls the sparsity level of the estimated condition specific networks.
}
  \item{lambda2}{
A nonnegative number. The hyperparameter controls the similarities between the estimated condition specific networks.
}
  \item{covType}{
A parameter to decide which approach we choose to compute the sample covariance matrices.

If covType = "pearson", it means that we compute (when input X represents data directly) the sample covariance matrices using Pearson correlation and that PNJGL is implemented based on Gaussian graphical models.  If data is normal, we suggest covType = "pearson".

If covType = "kendall" or "spearman", it means that we compute (when input X represents data directly) the sample covariance matrices using kendall's tau correlation or Spearman's rho correlation and that PNJGL is implemented based on nonparanormal graphical models. If data is non-normal, we suggest covType = "kendall" or "spearman".
}

\item{weights}{
Determines the putative sample size of each condition's data.  Allowed values: a vector with length two; "equal", giving each condition weight 1; "sample.size", giving each condition weight corresponding to its sample size.

}

\item{penalize.diagonal}{
Determines whether the sparsity penalty is applied to the diagonal of condition specific precision matrices.
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

\details{
Perturbed-node joint graphical lasso (PNJGL) is similar to FGL. The major difference between the two models is that FGL assumes that the differences between two condition specific networks arise from individual edges, while PNJGL assumes that network differences arise from certain nodes perturbed across conditions. A row-column overlap norm penalty is applied to the precision matrix difference to identify the crucial nodes driving network rewiring.

We solve the optimization problem of PNJGL using the ADMM algorithm. We accelerate the ADMM iterations by adaptively changing \eqn{\rho} in iterations. 
}

\value{
\item{Delta}{The estimateed difference of condition specific precision matrices.}
\item{Delta.graph.full}{The estimated differential network over all nodes.}
\item{Delta.graph.connected}{The estimated differential network over only the connected nodes.}
\item{Theta}{A list of the estimated condition specific precision matrices.}
\item{Theta.graph.full}{A list of the estimated condition specific networks over all nodes.}
\item{Theta.graph.connected}{A list of the estimated condition specific networks over only the connected nodes.}
}
\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Xiaohua Hu and Hong Yan (2017), DiffGraph: An R package for identifying gene
network rewiring using differential graphical models.\cr

Karthik Mohan, Palma London, Daniela Witten and Su-In Lee (2014). Node-based learning of multiple Gaussian graphicalmodels.
Journal of Machine Learning Research, 15(1), 445-488.
}

\author{
\packageAuthor{DiffGraph}

Maintainer: \packageMaintainer{DiffGraph}
}


\seealso{
\code{\link{Dtrace}}, \code{\link{FGL}}, \code{\link{pDNA}}, \code{\link{TCGA.BRCA}}, \code{\link{TCGA.GBM}}
}

\examples{

# Identify differential network between breast cancer subtypes
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
pnjgl.results= PNJGL(X, 0.3, 1.2, covType = "spearman")
net.pnjgl = pnjgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.pnjgl, vertex.size= degree(net.pnjgl)*.3, layout =layout_with_fr,
    vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange") 
# grab the coordinates from tkplot
l.pnjgl <- tkplot.getcoords(tkid) 
plot(net.pnjgl, layout=l.pnjgl,  vertex.size= degree(net.pnjgl)*.3,  vertex.color="red", 
    vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")



# Identify differential network between glioblastoma subtypes
data(TCGA.GBM)
X = TCGA.GBM$X[1,]
pnjgl.results= PNJGL(X, 0.2, 1.5, covType = "spearman")
net.pnjgl = pnjgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.pnjgl, vertex.size= degree(net.pnjgl)*.3, layout =layout_with_fr, 
    vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange") 
# grab the coordinates from tkplot
l.pnjgl <- tkplot.getcoords(tkid) 
plot(net.pnjgl, layout=l.pnjgl,  vertex.size= degree(net.pnjgl)*.3,  vertex.color="red", 
    vertex.label.cex =0.9,edge.width =1.5, edge.color="orange")
}

