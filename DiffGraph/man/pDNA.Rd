\name{pDNA}
\alias{pDNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior information-induced differential network analysis
}
\description{
This function is implemented to identify gene network rewiring using the pDNA model.
}
\usage{
pDNA(X, lambda, covType = "pearson", tol = 1e-05, maxiter = 500, rho = 0.1, 
    rho.incr = 1.05, rho.max = 1e+10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A list of input matrices. They can be data matrices or covariance matrices. If every matrix in the X is a symmetric matrix, the matrices are assumed to be covariance matrices among genes. If the input are data matrices, the rows represent the observations, and
   the columns represents the genes.
}
  \item{lambda}{
A nonnegative number. The hyperparameter controls the sparsity level of the estimated differential networks.
}
  \item{covType}{
A parameter to decide which approach we choose to compute the sample covariance matrices.

If covType = "pearson", it means that we compute (when input X represents data directly) the sample covariance matrices using Pearson correlation and that pDNA is implemented based on Gaussian graphical models.  If data is normal, we suggest covType = "pearson".

If covType = "kendall" or "spearman", it means that we compute (when input X represents data directly) the sample covariance matrices using kendall's tau correlation or Spearman's rho correlation and that pDNA is implemented based on nonparanormal graphical models. If data is non-normal, we suggest covType = "kendall" or "spearman".
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
pDNA is developed to deal with the sceneries where gene expression levels are measured using multiple data types. For each subject, we assume that the gene expression levels  are collected from \eqn{K} different data types. Let \eqn{S^{\left(kc\right)}} be the sample covariance matrix for the \eqn{K}th data types and the \eqn{g}th condition for \eqn{k=1,\ldots,K} and \eqn{g=1,2}, and let \eqn{\Delta^{\left(k\right)}} be the precision matrix difference for the \eqn{k}th data types. pDNA estimates \eqn{\Delta^{\left(1\right)}, \ldots, \Delta^{\left(K\right)}} by solving the following optimization problem:
\deqn{
 \begin{array}{ll}
\min\limits_{ \left\{ \Delta^k \right\}, \left\{ V^k \right\}} & \sum\limits_{k=1}^K  \left\{ \frac{1}{2} \left\langle   \Delta^2 ,  \frac{ S^{\left(k1\right)} S^{\left(k2\right)} + S^{\left(k2\right)}S^{\left(k1\right)} } {2 } \right\rangle - \left\langle \Delta, S^{\left(k1\right)} -  S^{\left(k2\right)} \right\rangle \right\}  + \lambda \sum\limits_{j=1}^p \sqrt{ \sum\limits_{i=1}^p \sqrt{\sum\limits_{k=1}^K \left| V^{\left(k\right)}_{ij} \right| } } \\[0.1in]
\mathrm{subject \ to} & \Delta^{\left(k\right)} = V^{\left(k\right)} + {\left( V^{\left(k\right)}\right)}^T,  \ k = 1, \ldots, K,
\end{array}
}
where the first term is the D-trace lasso used to infer precision matrix difference directly, the second term is a hierarchical group bridge penalty. The hierarchical group bridge penalty is developed to (i) integrate information across different data types, and (ii) capture important genes that drive the changes of network in a similar manner to PNJGL.\cr
}


\value{
\item{Delta}{List of the estimateed differences of condition specific precision matrices for different data types.}
\item{Delta.graph.full}{List of the estimated differential network over all nodes for different data types.}
\item{Delta.graph.connected}{List of the estimated differential network over only the connected nodes for different data types.}
\item{Delta.weight}{The difference of condition specific precision matrices averaged over different data types.}
\item{Delta.graph.weight.full}{The weighted differential network over all nodes, where the weight of a pair of genes is the frequency of data types in which they are identified as a differential edge.}
\item{Delta.graph.connected}{The weighted differential network over only the connected nodes, where the weight of a pair of genes is the frequency of data types in which they are identified as a differential edge.}
}
\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Xiaohua Hu and Hong Yan (2017), DiffGraph: An R package for identifying gene
network rewiring using differential graphical models.\cr

Xiao-Fei Zhang, Le Ou-Yang and Hong Yan (2017). Incorporating prior information into differential network
analysis using non-paranormal graphical models. Bioinformatics, 33(16), 2436-2445
}

\author{
\packageAuthor{DiffGraph}

Maintainer: \packageMaintainer{DiffGraph}
}


\seealso{
\code{\link{Dtrace}}, \code{\link{FGL}}, \code{\link{PNJGL}},  \code{\link{TCGA.BRCA}},\code{\link{TCGA.GBM}}
}


\examples{
  # 
  # # Identify differential network between breast cancer subtypes
  # data(TCGA.BRCA)
  # pdna.results= pDNA(TCGA.BRCA$X, 0.7, covType = "spearman")
  # net.pdna = pdna.results$Delta.graph.weight.connected
  # # Visualize the estimated differential network in an interactive manner. 
  # tkid <- tkplot(net.pdna, vertex.size= degree(net.pdna)*1.5, layout =layout_with_fr, 
  #                vertex.color="red", vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")
  # # Visualize the estimated differential network in a non-interactive manner.                
  # # grab the coordinates from tkplot
  # l.pdna <- tkplot.getcoords(tkid) 
  # plot(net.pdna, layout=l.pdna,  vertex.size= degree(net.pdna)*1.5,  vertex.color="red",
  #      vertex.label.cex=0.9,  edge.width =1.5, edge.color="orange")
  # 
  # 
  # # Identify differential network between glioblastoma subtypes
  # data(TCGA.GBM)
  # pdna.results= pDNA(TCGA.GBM$X, 0.5, covType = "spearman")
  # net.pdna = pdna.results$Delta.graph.weight.connected
  # # Visualize the estimated differential network in an interactive manner.
  # tkid <- tkplot(net.pdna, vertex.size= degree(net.pdna)*1.5, layout =layout_with_fr, 
  #                vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
  # # Visualize the estimated differential network in a non-interactive manner.
  # # grab the coordinates from tkplot
  # l.pdna <- tkplot.getcoords(tkid) 
  # plot(net.pdna, layout=l.pdna,  vertex.size= degree(net.pdna)*1.5,  vertex.color="red", 
  #      vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")
}

