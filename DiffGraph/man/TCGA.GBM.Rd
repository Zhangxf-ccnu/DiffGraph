\name{TCGA.GBM}
\alias{TCGA.GBM}
\docType{data}
\title{
TCGA glioblastoma data}
\description{
The TCGA glioblastoma gene expression data used in our case study. The data are obtained from the TCGA database. They are collected from three platforms: Agilent 244K Custom Gene Expression G450, Affymetrix HT Human Genome U133 Array Plate Set and Affymetrix Human Exon 1.0 ST Array.  They include expression measurments for 81 proneural tumors and 124 mesenchymal tumors. The data only include expression measurements of genes that overlap with the RTK/PI3K, p53, Rb signaling pathways which consists of the most frequently altered genes in glioblastoma.
}

\usage{data("TCGA.GBM")}

\source{
The Cancer Genome Atlas Research Network (2012), Comprehensive molecular portraits of human breast tumors. Nature. 490(7418), 61-70. (http://cancergenome.nih.gov/)
}
\references{
Xiao-Fei Zhang, Le Ou-Yang, Shuo Yang, Hong Qin, Xiaohua Hu and Hong Yan (2017), JGNI: A graphical model for joint gene network inference across multiple subpopulations and data types.\cr

Roel GW Verhaak, Katherine A Hoadley, Elizabeth Purdom, Victoria Wang, Yuan Qi, Matthew D Wilkerson, C Ryan
Miller, Li Ding, Todd Golub, Jill P Mesirov, et al. Integrated genomic analysis identifies clinically relevant subtypes of
glioblastoma characterized by abnormalities in pdgfra, idh1, egfr, and nf1. Cancer cell, 17(1):98-110, 2010.
}

\author{
\packageAuthor{DiffGraph}
}

\seealso{
  \code{\link{TCGA.BRCA}}
}


\examples{
data(TCGA.GBM)
## maybe str(TCGA.GBM) ...
}
