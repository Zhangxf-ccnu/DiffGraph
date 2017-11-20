library("DiffGraph")

rm(list = ls())

# Identify differential network between breast cancer subtypes using Dtrace with covType = "spearman"

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



# Identify differential network between breast cancer subtypes using Dtrace with covType = "kendall"
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
dtrace.results= Dtrace(X, 0.45, covType = "kendall")
net.dtrace = dtrace.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.dtrace, vertex.size= degree(net.dtrace)*1.5, layout =layout_with_fr,
               vertex.color="red", vertex.label.cex=0.8, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.dtrace <- tkplot.getcoords(tkid)
plot(net.dtrace, layout=l.dtrace,  vertex.size= degree(net.dtrace)*1.5,  vertex.color="red", 
     vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")




# Identify differential network between breast cancer subtypes using Dtrace with covType = "pearson"
data(TCGA.GBM)
X = TCGA.GBM$X[1,]
dtrace.results= Dtrace(X, 0.45, covType = "pearson")
net.dtrace = dtrace.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.dtrace, vertex.size= degree(net.dtrace)*1.5, layout =layout_with_fr,
               vertex.color="red", vertex.label.cex=0.8, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.dtrace <- tkplot.getcoords(tkid)
plot(net.dtrace, layout=l.dtrace,  vertex.size= degree(net.dtrace)*1.5,  vertex.color="red", 
     vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")




# Identify differential network between glioblastoma subtypes using Dtrace with covType = "spearman"
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





##Identify differential network between breast cancer subtypes using FGL with covType = "spearman"
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
fgl.results= FGL(X, 0.05, 0.2, covType = "spearman") 
net.fgl = fgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner
tkid <- tkplot(net.fgl, vertex.size= degree(net.fgl)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
# grab the coordinates from tkplot
l.fgl <- tkplot.getcoords(tkid) 
# Visualize the estimated differential network in a non-interactive manner.
plot(net.fgl, layout=l.fgl, vertex.size= degree(net.fgl)*1.5,  vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")


##Identify differential network between breast cancer subtypes using FGL with covType = "kendall"
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
fgl.results= FGL(X, 0.05, 0.2, covType = "kendall") 
net.fgl = fgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner
tkid <- tkplot(net.fgl, vertex.size= degree(net.fgl)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
# grab the coordinates from tkplot
l.fgl <- tkplot.getcoords(tkid) 
# Visualize the estimated differential network in a non-interactive manner.
plot(net.fgl, layout=l.fgl, vertex.size= degree(net.fgl)*1.5,  vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")


##Identify differential network between breast cancer subtypes using FGL with covType = "pearson"
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
fgl.results= FGL(X, 0.05, 0.2, covType = "pearson") 
net.fgl = fgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner
tkid <- tkplot(net.fgl, vertex.size= degree(net.fgl)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
# grab the coordinates from tkplot
l.fgl <- tkplot.getcoords(tkid) 
# Visualize the estimated differential network in a non-interactive manner.
plot(net.fgl, layout=l.fgl, vertex.size= degree(net.fgl)*1.5,  vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")




## Identify differential network between glioblastoma subtypes using FGL with covType = "spearman"
data(TCGA.GBM)
X = TCGA.GBM$X[1,]
fgl.results= FGL(X, 0.1, 0.16, covType = "spearman")
net.fgl = fgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner
tkid <- tkplot(net.fgl, vertex.size= degree(net.fgl)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange")
# grab the coordinates from tkplot
l.fgl <- tkplot.getcoords(tkid) 
# Visualize the estimated differential network in a non-interactive manner.
plot(net.fgl, layout=l.fgl, vertex.size= degree(net.fgl)*1.5, vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")




# Identify differential network between breast cancer subtypes using PNJGL with covType = "spearman"
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




# Identify differential network between breast cancer subtypes using PNJGL with covType = "kendall"
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
pnjgl.results= PNJGL(X, 0.3, 1.2, covType = "kendall")
net.pnjgl = pnjgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.pnjgl, vertex.size= degree(net.pnjgl)*.3, layout =layout_with_fr,
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange") 
# grab the coordinates from tkplot
l.pnjgl <- tkplot.getcoords(tkid) 
plot(net.pnjgl, layout=l.pnjgl,  vertex.size= degree(net.pnjgl)*.3,  vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")




# Identify differential network between breast cancer subtypes using PNJGL with covType = "pearson"
data(TCGA.BRCA)
X = TCGA.BRCA$X[1,]
pnjgl.results= PNJGL(X, 0.3, 1.2, covType = "pearson")
net.pnjgl = pnjgl.results$Delta.graph.connected
# Visualize the estimated differential network in an interactive manner.
tkid <- tkplot(net.pnjgl, vertex.size= degree(net.pnjgl)*.3, layout =layout_with_fr,
               vertex.color="red", vertex.label.cex=0.9, edge.width =1.5, edge.color="orange") 
# grab the coordinates from tkplot
l.pnjgl <- tkplot.getcoords(tkid) 
plot(net.pnjgl, layout=l.pnjgl,  vertex.size= degree(net.pnjgl)*.3,  vertex.color="red", 
     vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")





# Identify differential network between glioblastoma subtypes using PNJGL with covType = "spearman"
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





# Identify differential network between breast cancer subtypes using pDNA with covType = "spearman"
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




# Identify differential network between breast cancer subtypes using pDNA with covType = "kendall"
data(TCGA.BRCA)
pdna.results= pDNA(TCGA.BRCA$X, 0.9, covType = "kendall")
net.pdna = pdna.results$Delta.graph.weight.connected
# Visualize the estimated differential network in an interactive manner. 
tkid <- tkplot(net.pdna, vertex.size= degree(net.pdna)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.pdna <- tkplot.getcoords(tkid) 
plot(net.pdna, layout=l.pdna,  vertex.size= degree(net.pdna)*1.5,  vertex.color="red",
     vertex.label.cex=0.9,  edge.width =1.5, edge.color="orange")





# Identify differential network between breast cancer subtypes using pDNA with covType = "pearson"
data(TCGA.BRCA)
pdna.results= pDNA(TCGA.BRCA$X, 1.1, covType = "pearson")
net.pdna = pdna.results$Delta.graph.weight.connected
# Visualize the estimated differential network in an interactive manner. 
tkid <- tkplot(net.pdna, vertex.size= degree(net.pdna)*1.5, layout =layout_with_fr, 
               vertex.color="red", vertex.label.cex =0.9, edge.width =1.5, edge.color="orange")
# Visualize the estimated differential network in a non-interactive manner.                
# grab the coordinates from tkplot
l.pdna <- tkplot.getcoords(tkid) 
plot(net.pdna, layout=l.pdna,  vertex.size= degree(net.pdna)*1.5,  vertex.color="red",
     vertex.label.cex=0.9,  edge.width =1.5, edge.color="orange")




# Identify differential network between glioblastoma subtypes using pDNA with covType = "spearman"
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
