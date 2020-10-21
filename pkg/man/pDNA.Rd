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
pDNA(X, lambda, covType = "pearson", tol = 1e-05, maxiter = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A matrix (K x 2) of list of input matrices, where K is the number of data types. The (k, g)-th element is the data matrix 
  corresponding to the k-th data type and the g-th condition. They can be data matrices or covariance matrices. If every matrix in   X is a symmetric matrix, the matrices are assumed to be covariance matrices among genes. If the input are data matrices, the rows represent the samples, and the columns represents the genes.
}
  \item{lambda}{
A nonnegative number. The tuning controls the sparsity level of the estimated differential networks.
}
  \item{covType}{
A parameter to decide which approach we choose to compute the sample covariance matrices.

If covType = "pearson", it means that we compute (when input X represents data) the sample covariance matrices using Pearson correlation and that pDNA is implemented based on Gaussian graphical models.  If data is normal, we suggest covType = "pearson".

If covType = "kendall" or "spearman", it means that we compute (when input X represents data) the sample covariance matrices using kendall's tau correlation or Spearman's rho correlation and that pDNA is implemented based on nonparanormal graphical models. If data is non-normal, we suggest covType = "kendall" or "spearman".
}
  \item{tol}{
The tolerance parameter for convergence criteria.
}
  \item{maxiter}{
The maximum number of iterations.
}
}
\details{
Prior information-induced differential network analysis (pDNA) is an extension of Dtrace to deal with the settings where gene expression measurements are collected using multiple data types. Information across different data types are integrated. Furthermore, pDNA can capture important genes that drive the changes of network in a similar manner to PNJGL.

We solve the optimization problem of pDNA using the local linear approximation method and the accelerated proximal gradient method. 
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
 
# Identify differential network between breast cancer subtypes
data(TCGA.BRCA)
pdna.results= pDNA(TCGA.BRCA$X, 0.9, covType = "spearman")
net.pdna = pdna.results$Delta.graph.weight.connected
# Visualize the estimated differential network in an interactive manner. 
tkid <- tkplot(net.pdna, vertex.size= degree(net.pdna)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.pdna <- tkplot.getcoords(tkid) 
plot(net.pdna, layout=l.pdna,  vertex.size= degree(net.pdna)*1.5,  vertex.color="red",
     vertex.label.cex=0.9,  edge.width =1.5, edge.color="orange")


# Identify differential network between glioblastoma subtypes
data(TCGA.GBM)
pdna.results= pDNA(TCGA.GBM$X, 0.64, covType = "spearman")
net.pdna = pdna.results$Delta.graph.weight.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.pdna, vertex.size= degree(net.pdna)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.
# grab the coordinates from tkplot
l.pdna <- tkplot.getcoords(tkid) 
plot(net.pdna, layout=l.pdna,  vertex.size= degree(net.pdna)*1.5,  vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")
}
