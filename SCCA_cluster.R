####### hierarchical clustering -----
####### clear the environment
rm (list = ls(all.names=TRUE))
####### install packages
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
install.packages("dendextend")
install.packages("https://cran.rstudio.com/bin/windows/contrib/4.2/rlang_1.1.0.zip",repo=NULL, type="source")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("viridis")
install.packages("circlize")             
install.packages("RColorBrewer")        
install.packages("factoextra")
install.packages("cluster")
install.packages("rstatix")

library(R.matlab)
library(rlang)
library(dplyr)
library(ComplexHeatmap)
library(dendextend)
library(ggplot2)
library(viridis)
library(circlize)
library(RColorBrewer)   
library(tidyverse)
library(cluster)
library(factoextra)
library(rstatix)
library(stats)


rm(list=ls())
#######  clustering for results produced by R
outputPath<-".../output/"
# load sCCA canonical variate
setwd(".../sCCAresults/")
load("scca_cors.Rdata")
load("scca_perm10000.Rdata")
Xu.use<-scca.cand$Xu[,1:2]    # here is the canonical variate for brain imaging that pass the permutation test (< FDR 0.05)
data<-as.matrix(Xu.use)

#### ---- Settings ---- ####
set.seed(874) 

## do hierarchical clustering analysis
# basic info
str(data)
summary(data)
any(is.na(data))

data.scale<-as.data.frame(scale(data))
summary(data.scale)


#### ---- Evaluate linkage methods ---- ####
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
ac <- function(x) {
  agnes(data.scale, method = x)$ac
}                 # chapter 5 of Kaufman and Rousseeuw (1990)
map_dbl(m, ac)   

                  
### ---- Run hierarchical clustering ---- ####
dist.mat <- dist(data.scale, method = 'euclidean')
hclust <- hclust(dist.mat, method = 'ward.D')       
p.den<-plot(hclust)
# cut the dendrogram in order to create the desired number of clusters
k=3
cut <- cutree(hclust, k = k)
dend.obj <- as.dendrogram(hclust)
col.dend <- dendextend::color_branches(dend.obj, h = k)
p.den<-plot(col.dend)


#### ---- Evaluate optimal clustering solution ---- ####
# determine optimal n clusters : dendrogram
pdf(paste0(outputPath,"cluster_den.pdf"), paper='a4')
p.den<-plot(hclust, hang = -1, cex = 0.5)
p.den
dev.off()

# determine optimal n clusters : elbow method
pdf(paste0(outputPath,"cluster_TotalWithinSumSquare.pdf"), paper='a4')
p.sumSqu<-fviz_nbclust(data.scale, FUN = hcut, method = "wss")
p.sumSqu
dev.off()

# determine optimal n clusters : average silhouette method
pdf(paste0(outputPath,"cluster_avgSilhouette.pdf"), paper='a4')
p.Sil<-fviz_nbclust(data.scale, FUN = hcut, method = "silhouette")
p.Sil
dev.off()


#### ---- Select clusters and save data ---- ####
cut <- cutree(hclust, k = 3)
save(list = c("hclust","cut"), file = paste0(outputPath,"cluster_results_Cluster3.Rdata"))

Cidx<-cut
Cidx<-as.data.frame(Cidx)
names(Cidx)<-"Cidx"
write.csv(Cidx,file=paste0(outputPath,"Cidx_cluster3.csv"),row.names=FALSE)

table(Cidx['Cidx'])


#### ---- plot the heatmap with dendrogram ---- ####
mycol <- colorRamp2(c(-3, 0, 3),c("#1f70a9","#EEEEEE","#ea7827"))
row_dend = hclust(dist(as.matrix(data.scale),method = "euclidean"), method = "ward.D")  
grow_dend = color_branches(row_dend, k = 3)

pdf(paste0(outputPath,"cluster_heatmap_cluster3.pdf"), paper='a4')
p.heatmap<-Heatmap(as.matrix(data.scale), 
                   col = mycol,
                   color_space = "LAB",
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "ward.D",         
                   cluster_columns = FALSE,
                   heatmap_height = unit(10, "cm"),
                   heatmap_width = unit(7, "cm"),
                   row_dend_width = unit(4, "cm"),
                   split = 3,
                   cluster_rows = row_dend,
                   row_dend_gp=gpar(lwd = 1.2))  
p.heatmap
dev.off()


