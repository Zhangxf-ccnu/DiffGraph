\name{TCGA.BRCA}
\alias{TCGA.BRCA}
\docType{data}
\title{
TCGA breast cancer data}
\description{
The TCGA breast cancer gene expression datasets used in our case study. The datasets are obtained from the TCGA database. They are collected from three platforms: mRNA expression (Agilent G450 microarray), mRNA expression (RNA sequencing) and copy number variants (Affymetrix genome-wide human SNP Array 6.0). They include expression profiles for 207 luminal A cancers and 86 basal-like cancers.
The data only include expression measurement of genes that overlap with the breast cancer pathway collected from the Kyoto Encyclopedia of Genes and Genomes database.
}

\usage{data("TCGA.BRCA")}


\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Xiaohua Hu and Hong Yan (2017), DiffGraph: An R package for identifying gene
network rewiring using differential graphical models. \cr

The Cancer Genome Atlas Research Network (2012), Comprehensive molecular portraits of human breast tumors. Nature. 490(7418), 61-70. (http://cancergenome.nih.gov/)
}


\author{
\packageAuthor{DiffGraph}
}

\seealso{
\code{\link{Dtrace}}, \code{\link{FGL}}, \code{\link{PNJGL}}, \code{\link{pDNA}},  \code{\link{TCGA.GBM}}
}


\examples{
data(TCGA.BRCA)
## maybe str(TCGA.BRCA) ...
}

