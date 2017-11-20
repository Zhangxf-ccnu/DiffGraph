README file for R package supporting the paper "DiffGraph: An R package for identifying gene network rewiring using differential graphical models".


Contents of this archive
------------------------
This archive contains 
(1) DiffGraph_1.1.1.tar: Package source. 
(2) DiffGraph-manual.pdf: Reference manual.
(3) example.R: Examples for step by step usages for the DiffGraph package.

The DiffGraph pacakge depends on the following existing packages. To use them, simply
library('igraph')
library('MASS')
library('Matrix')
If you don't have these packages, simply use
install.packages("igraph")
install.packages("MASS")
install.packages("Matrix")

Install the DiffGraph pacakge:
Download the package source from https://github.com/Zhangxf-ccnu/DiffGraph 
and then install the package using the following command:
install.packages("path/DiffGraph_1.1.1.tar.gz", type = "source")
where "path" is the path where the "DiffGraph_1.1.1.tar.gz" is located.

Useage
Load the library DiffGraph in R console, by running:
library("DiffGraph")

Simply run the one of diffential graphical models on your favorite datasets. For example,
data(TCGA.BRCA)
pdna.results= pDNA(TCGA.BRCA$X, 0.9, covType = "spearman")

For detialed usages, please refer to "DiffGraph-manual.pdf".
For more examples, please refer to "example.R"

Please do not hesitate to contact Dr. Xiao-Fei Zhang at zhangxf@mail.ccnu.edu.cn to 
seek any clarifications regarding any contents or operation of the archive.
