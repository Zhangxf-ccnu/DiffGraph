README file for R package supporting the paper: Xiao-Fei Zhang, Le Ou-Yang*, Shuo Yang, Xiaohua Hu, Hong Yan, DiffGraph: an R package for identifying gene network rewiring using differential graphical models, Bioinformatics, 2018,34(9): 1571-1573. 

Contents of this archive
------------------------
This archive contains 
(1) pkg: subdirectory that contains the R package. 
(2) DiffGraph-manual.pdf: Reference manual.
(3) example.R: Examples for step by step usages for the DiffGraph package.

The DiffGraph package has the following R-package dependencies: igraph, MASS, Matrix. The dependent packages will be automatically installed 
along with DiffGraph. You can use the following commands to install DiffGraph from GitHub.

Step 1. If the devtools package has been not installed, install the devtools package first. Invoke R and then type
install.packages("devtools")

Step 2. Load the devtools package.
library("devtools")

Step 3. Install the DiffGraph package from GitHub.
install_github("Zhangxf-ccnu/DiffGraph", subdir="pkg")


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
