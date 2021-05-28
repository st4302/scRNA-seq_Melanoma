library(GSVA)
library(CHETAH)
library(SingleR)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(Seurat)
library(msigdbr)
library(readr)
library(celldex)
library(SingleR)
library(scuttle)
library(scran)
library(nichenetr)
library(tidyverse)

memory.limit()
memory.limit(size=65000)

###SINGLER####


########Resolution 1
#Current
dat <- readRDS("seurat_obj.rds")
dim(dat) # 2000x18474
# Get cell identity classes
Idents(object = dat)
levels(x = dat) #C01, C02 ... C93
# View metadata data frame, stored in object@meta.data
dat[[]]


Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

#Annotate each cell via SingleR() function. Assignment scores based on Spearman correlation across markers
pred.hesc <- SingleR(test = dat@assays$RNA@data, ref = Mouseseq, assay.type.test=1, labels = Mouseseq$label.main)

head(pred.hesc) #DataFrame with 18474 rows and 5 columns
colnames(pred.hesc)
#[1] "scores"        "first.labels"  "tuning.scores" "labels"       
#[5] "pruned.labels"


dat[["SingleR.labels"]] <- pred.hesc$labels
Idents(dat) <- "SingleR.labels"
png("p25.png", height = 5 , width = 7, units="in", res=150)
DimPlot(dat, reduction = "umap")


dat[["gsva_res.gene.symbol"]] <- gsva_res$gene.symbol
Idents(dat) <- "SingleR.labels"
png("p25.png", height = 5 , width = 7, units="in", res=150)
DimPlot(dat, reduction = "umap")

table(pred.hesc$labels)

#Adipocytes        Astrocytes           B cells    Cardiomyocytes
#1041                 2                14                 6
#Endothelial cells  Epithelial cells      Erythrocytes       Fibroblasts
#1757               445               118             13956
#Granulocytes       Hepatocytes       Macrophages         Microglia
#1               498               250                 1
#Monocytes           Neurons          NK cells  Oligodendrocytes
#106                88                26               102
#T cells
#63



png("p1.png", height = 5 , width = 7, units="in", res=150)
p1 =plotScoreHeatmap(pred.hesc)


png("p2.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.hesc, ncol=9, size=0.5) #many warnings but plot
#Warning messages:
#  1: In max(data$density) : no non-missing arguments to max; returning -Inf
#2: Computation failed in `stat_ydensity()`:
#  replacement has 1 row, data has 0
#3: In max(data$density) : no non-missing arguments to max; returning -Inf
#4: Computation failed in `stat_ydensity()`:
#  replacement has 1 row, data has 0
#5: In max(data$density) : no non-missing arguments to max; returning -Inf
#6: Computation failed in `stat_ydensity()`:
#  replacement has 1 row, data has 0

summary(is.na(pred.hesc$pruned.labels))

#Mode   FALSE    TRUE
#logical   18473       1

#Aith a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 


pred.grun <- SingleR(test=dat@assays$RNA@data, ref=Mouseseq, labels=Mouseseq$label.main, de.method="wilcox") 


png("p3.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.grun)

png("p4.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.grun, ncol = 9, size=0.5) #no warnings

table(pred.grun$labels)

Adipocytes        Astrocytes           B cells    Cardiomyocytes
1322                40               340               261
Dendritic cells Endothelial cells  Epithelial cells      Erythrocytes
256              2064              1182               240
Fibroblasts      Granulocytes       Hepatocytes       Macrophages
8232               408              1531              1073
Microglia         Monocytes           Neurons          NK cells
76               149                74                98
Oligodendrocytes           T cells
111              1017



summary(is.na(pred.grun$pruned.labels)) 

#Mode   FALSE    TRUE
#logical   18472       2

########Resolution 1 Labeld fine
#Current
dat <- readRDS("seurat_obj.rds")
dim(dat) # 2000x18474

Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

#Annotate each cell via SingleR() function. Assignment scores based on Spearman correlation across markers
pred.hesc2 <- SingleR(test = dat@assays$RNA@data, ref = Mouseseq, assay.type.test=1, labels = Mouseseq$label.fine)

dat[["SingleR1.labels"]] <- pred.hesc2$labels
Idents(dat) <- "SingleR1.labels"
png("p25.png", height = 5 , width = 7, units="in", res=150)
DimPlot(dat, reduction = "umap")

pred.hesc2 #DataFrame with 18474 rows and 5 columns

table(pred.hesc2$labels)

#Adipocytes                 aNSCs  Astrocytes activated
#497                    70                     2
#B cells        Cardiomyocytes     Endothelial cells
#6                     7                  2034
#Ependymal          Erythrocytes           Fibroblasts
#331                    70                 10707
#Fibroblasts activated Fibroblasts senescent          Granulocytes
#892                  1494                     1
#Hepatocytes           Macrophages Macrophages activated
#736                    15                   126
#Microglia   Microglia activated             Monocytes
#1                    26                    45
#NK cells                  NPCs      Oligodendrocytes
#34                  1223                    81
#OPCs                 qNSCs               T cells
#20                    39                    17



png("p5.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.hesc2)  #warnings

png("p6.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.hesc2, ncol=9, size=0.5)

summary(is.na(pred.hesc2$pruned.labels))

#Mode   FALSE    TRUE
#logical   18472       2

#With a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 


pred.grun2 <- SingleR(test=dat@assays$RNA@data, ref=Mouseseq, labels=Mouseseq$label.fine, de.method="wilcox") 


png("p7.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.grun2)

png("p8.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.grun2, ncol = 9, size=0.5)

table(pred.grun2$labels)

#Adipocytes                 aNSCs            Astrocytes
#1379                   209                    25
#Astrocytes activated               B cells        Cardiomyocytes
#6                   402                   460
#Dendritic cells     Endothelial cells             Ependymal
#167                   894                   957
#Erythrocytes           Fibroblasts Fibroblasts activated
#277                  6317                   160
#Fibroblasts senescent          Granulocytes           Hepatocytes
#2546                   499                  1746
#Macrophages Macrophages activated             Microglia
#383                   338                    11
#Microglia activated             Monocytes               Neurons
#52                   155                     3
#NK cells                  NPCs      Oligodendrocytes
#66                   264                   171
#OPCs                 qNSCs               T cells
#97                   154                   736

summary(is.na(pred.grun2$pruned.labels)) 

#Mode   FALSE    TRUE
#logical   18467      7




########Resolution 1 ##prob correct
#Current
dat <- readRDS("seurat_obj.rds")
dim(dat) # 2000x18474
Idents(dat) <- "snn_res.1"
#Idents(dat) <- "SingleR.labels"
# View metadata data frame, stored in object@meta.data
dat[[]] #not sure if it chooses it



Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

#Annotate each cell via SingleR() function. Assignment scores based on Spearman correlation across markers
pred.hesc <- SingleR(test = dat@assays$RNA@data, ref = Mouseseq, assay.type.test=1, labels = Mouseseq$label.main)

#dat <- AddMetaData(object = dat, metadata = labels)
dat[["SingleR1.labels"]] <- pred.hesc$labels
Idents(dat) <- "SingleR1.labels"
png("p25.png", height = 5 , width = 7, units="in", res=150)
DimPlot(dat, reduction = "umap")
DimPlot(dat, reduction = "umap", group.by="snn_res.1")
dat[[]]
gsva_res[[]]


pred.hesc #DataFrame with 18474 rows and 5 columns

table(pred.hesc$labels)

#Adipocytes        Astrocytes           B cells    Cardiomyocytes
#1041                 2                14                 6
#Endothelial cells  Epithelial cells      Erythrocytes       Fibroblasts
#1757               445               118             13956
#Granulocytes       Hepatocytes       Macrophages         Microglia
#1               498               250                 1
#Monocytes           Neurons          NK cells  Oligodendrocytes
#106                88                26               102
#T cells
#63



png("p9.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.hesc)

png("p10.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.hesc, ncol=9, size=0.5) #many warnings but plot


summary(is.na(pred.hesc$pruned.labels))

#Mode   FALSE    TRUE
#logical   18473       1

#Aith a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 


pred.grun <- SingleR(test=dat@assays$RNA@data, ref=Mouseseq, labels=Mouseseq$label.main, de.method="wilcox") 


png("p11.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.grun)

png("p12.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.grun, ncol = 9, size=0.5) #no warnings

table(pred.grun$labels)

Adipocytes        Astrocytes           B cells    Cardiomyocytes
1322                40               340               261
Dendritic cells Endothelial cells  Epithelial cells      Erythrocytes
256              2064              1182               240
Fibroblasts      Granulocytes       Hepatocytes       Macrophages
8232               408              1531              1073
Microglia         Monocytes           Neurons          NK cells
76               149                74                98
Oligodendrocytes           T cells
111              1017



summary(is.na(pred.grun$pruned.labels)) 

#Mode   FALSE    TRUE
#logical   18472       2



########Resolution 1 Labeld fine
#Current
dat <- readRDS("seurat_obj.rds")
dim(dat) # 2000x18474

Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

#Annotate each cell via SingleR() function. Assignment scores based on Spearman correlation across markers
pred.hesc2 <- SingleR(test = dat@assays$RNA@data, ref = Mouseseq, assay.type.test=1, labels = Mouseseq$label.fine)

pred.hesc2 #DataFrame with 18474 rows and 5 columns

table(pred.hesc2$labels)

dat[["SingleR1.labels"]] <- pred.hesc$labels
Idents(dat) <- "SingleR.labels"
png("p25.png", height = 5 , width = 7, units="in", res=150)
DimPlot(dat, reduction = "umap")

#Adipocytes                 aNSCs  Astrocytes activated
#497                    70                     2
#B cells        Cardiomyocytes     Endothelial cells
#6                     7                  2034
#Ependymal          Erythrocytes           Fibroblasts
#331                    70                 10707
#Fibroblasts activated Fibroblasts senescent          Granulocytes
#892                  1494                     1
#Hepatocytes           Macrophages Macrophages activated
#736                    15                   126
#Microglia   Microglia activated             Monocytes
#1                    26                    45
#NK cells                  NPCs      Oligodendrocytes
#34                  1223                    81
#OPCs                 qNSCs               T cells
#20                    39                    17



png("p13.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.hesc2)  #warnings

png("p14.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.hesc2, ncol=9, size=0.5)

summary(is.na(pred.hesc2$pruned.labels))

#Mode   FALSE    TRUE
#logical   18472       2

#With a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 


pred.grun2 <- SingleR(test=dat@assays$RNA@data, ref=Mouseseq, labels=Mouseseq$label.fine, de.method="wilcox") 


png("p15.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.grun2)

png("p16.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.grun2, ncol = 9, size=0.5)

table(pred.grun2$labels)

#Adipocytes                 aNSCs            Astrocytes
#1379                   209                    25
#Astrocytes activated               B cells        Cardiomyocytes
#6                   402                   460
#Dendritic cells     Endothelial cells             Ependymal
#167                   894                   957
#Erythrocytes           Fibroblasts Fibroblasts activated
#277                  6317                   160
#Fibroblasts senescent          Granulocytes           Hepatocytes
#2546                   499                  1746
#Macrophages Macrophages activated             Microglia
#383                   338                    11
#Microglia activated             Monocytes               Neurons
#52                   155                     3
#NK cells                  NPCs      Oligodendrocytes
#66                   264                   171
#OPCs                 qNSCs               T cells
#97                   154                   736

summary(is.na(pred.grun2$pruned.labels)) 

#Mode   FALSE    TRUE
#logical   18467      7


########Resolution 0.4 ##prob correct
#Current
dat <- readRDS("seurat_obj.rds")
dim(dat) # 2000x18474
Idents(dat) <- "snn_res.0.4"


Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

#Annotate each cell via SingleR() function. Assignment scores based on Spearman correlation across markers
pred.hesc <- SingleR(test = dat@assays$RNA@data, ref = Mouseseq, assay.type.test=1, labels = Mouseseq$label.main)

pred.hesc #DataFrame with 18474 rows and 5 columns

table(pred.hesc$labels)

#Adipocytes        Astrocytes           B cells    Cardiomyocytes
#1041                 2                14                 6
#Endothelial cells  Epithelial cells      Erythrocytes       Fibroblasts
#1757               445               118             13956
#Granulocytes       Hepatocytes       Macrophages         Microglia
#1               498               250                 1
#Monocytes           Neurons          NK cells  Oligodendrocytes
#106                88                26               102
#T cells
#63



png("p17.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.hesc)

png("p18.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.hesc, ncol=9, size=0.5) #many warnings but plot
Warning messages:
  1: In max(data$density) : no non-missing arguments to max; returning -Inf
2: Computation failed in `stat_ydensity()`:
  replacement has 1 row, data has 0
3: In max(data$density) : no non-missing arguments to max; returning -Inf
4: Computation failed in `stat_ydensity()`:
  replacement has 1 row, data has 0
5: In max(data$density) : no non-missing arguments to max; returning -Inf
6: Computation failed in `stat_ydensity()`:
  replacement has 1 row, data has 0

summary(is.na(pred.hesc$pruned.labels))

#Mode   FALSE    TRUE
#logical   18473       1

#Aith a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 

#not done

pred.grun <- SingleR(test=dat@assays$RNA@data, ref=Mouseseq, labels=Mouseseq$label.main, de.method="wilcox") 


png("p19.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.grun)

png("p20.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.grun, ncol = 9, size=0.5) #no warnings

table(pred.grun$labels)

Adipocytes        Astrocytes           B cells    Cardiomyocytes
1322                40               340               261
Dendritic cells Endothelial cells  Epithelial cells      Erythrocytes
256              2064              1182               240
Fibroblasts      Granulocytes       Hepatocytes       Macrophages
8232               408              1531              1073
Microglia         Monocytes           Neurons          NK cells
76               149                74                98
Oligodendrocytes           T cells
111              1017



summary(is.na(pred.grun$pruned.labels)) 

#Mode   FALSE    TRUE
#logical   18472       2

########Resolution 1 Labeld fine
#Current
dat <- readRDS("seurat_obj.rds")
dim(dat) # 2000x18474

Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

#Annotate each cell via SingleR() function. Assignment scores based on Spearman correlation across markers
pred.hesc2 <- SingleR(test = dat@assays$RNA@data, ref = Mouseseq, assay.type.test=1, labels = Mouseseq$label.fine)

pred.hesc2 #DataFrame with 18474 rows and 5 columns

table(pred.hesc2$labels)

#Adipocytes                 aNSCs  Astrocytes activated
#497                    70                     2
#B cells        Cardiomyocytes     Endothelial cells
#6                     7                  2034
#Ependymal          Erythrocytes           Fibroblasts
#331                    70                 10707
#Fibroblasts activated Fibroblasts senescent          Granulocytes
#892                  1494                     1
#Hepatocytes           Macrophages Macrophages activated
#736                    15                   126
#Microglia   Microglia activated             Monocytes
#1                    26                    45
#NK cells                  NPCs      Oligodendrocytes
#34                  1223                    81
#OPCs                 qNSCs               T cells
#20                    39                    17



png("p21.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.hesc2)  #warnings

png("p22.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.hesc2, ncol=9, size=0.5)

summary(is.na(pred.hesc2$pruned.labels))

#Mode   FALSE    TRUE
#logical   18472       2

#With a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 


pred.grun2 <- SingleR(test=dat@assays$RNA@data, ref=Mouseseq, labels=Mouseseq$label.fine, de.method="wilcox") 


png("p23.png", height = 5 , width = 7, units="in", res=150)
plotScoreHeatmap(pred.grun2)

png("p24.png", height = 15 , width = 17, units="in", res=150)
plotDeltaDistribution(pred.grun2, ncol = 9, size=0.5)

table(pred.grun2$labels)

#Adipocytes                 aNSCs            Astrocytes
#1379                   209                    25
#Astrocytes activated               B cells        Cardiomyocytes
#6                   402                   460
#Dendritic cells     Endothelial cells             Ependymal
#167                   894                   957
#Erythrocytes           Fibroblasts Fibroblasts activated
#277                  6317                   160
#Fibroblasts senescent          Granulocytes           Hepatocytes
#2546                   499                  1746
#Macrophages Macrophages activated             Microglia
#383                   338                    11
#Microglia activated             Monocytes               Neurons
#52                   155                     3
#NK cells                  NPCs      Oligodendrocytes
#66                   264                   171
#OPCs                 qNSCs               T cells
#97                   154                   736

summary(is.na(pred.grun2$pruned.labels)) 

#Mode   FALSE    TRUE
#logical   18467      7
