###ssGSEA


suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
library(ComplexHeatmap)

dat <- readRDS("seurat_obj.rds")
ES <- readRDS("ES.rds")
sce <- readRDS("sce.rds")

#GS <- getGeneSets(species = "Mus musculus", library = "C4")
GS <- getGeneSets(species = "Mus musculus", library = "C7")

sce <- as.SingleCellExperiment(dat)
counts <-GetAssayData(object=dat)
umap <- reducedDim(sce,"UMAP")
input_mel <- SingleCellExperiment(assays = list(counts = counts), reducedDims= SimpleList(UMAP=umap))
ES <- enrichIt(obj = input_mel, gene.sets = GS, groups = 1000, cores = 8)
#ES.sce <- enrichIt(obj = pbmc.sce, gene.sets = GS, groups = 1000, cores = 2)

#write.table(ES, file = "ES.txt", sep = "\t",
#            row.names = TRUE, col.names = NA)
#saveRDS(ES, file = "ES.rds")


dat <- AddMetaData(dat, ES)
dat@meta.data$active.idents <- dat@active.ident








ES2 <- data.frame(dat[[]], Idents(dat))
colnames(ES2)[ncol(ES2)] <- "cluster"
write.table(ES2, file = "ES2.tsv", sep = "\t", quote = FALSE)

#PCA
PCA <- performPCA(enriched = ES2, groups = c("Type", "cluster"))
png("pca.png", height = 15 , width = 17, units="in", res=150)
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)

#Top 10 contributors
png("pcaplot.png", height = 15 , width = 17, units="in", res=150)
masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))


#Heatmap Top 10 contributors
dittoHeatmap(dat, genes = NULL, 
             metas = c("KAECH_DAY8_EFF_VS_DAY15_EFF_CD8_TCELL_UP", "GSE36476_CTRL_VS_TSST_ACT_40H_MEMORY_CD4_TCELL_OLD_DN", 
                       "GSE14415_TCONV_VS_FOXP3_KO_INDUCED_TREG_DN", "GSE24634_TEFF_VS_TCONV_DAY5_IN_CULTURE_UP",
                       "GSE36476_CTRL_VS_TSST_ACT_72H_MEMORY_CD4_TCELL_YOUNG_DN", "GSE36476_CTRL_VS_TSST_ACT_40H_MEMORY_CD4_TCELL_YOUNG_DN",
                       "GSE37532_WT_VS_PPARG_KO_VISCERAL_ADIPOSE_TISSUE_TREG_UP", "GSE30962_PRIMARY_VS_SECONDARY_ACUTE_LCMV_INF_CD8_TCELL_UP",
                       "GSE13547_CTRL_VS_ANTI_IGM_STIM_BCELL_2H_UP", "GSE24634_TEFF_VS_TCONV_DAY7_IN_CULTURE_UP"), 
             heatmap.colors = rev(colorblind_vector(50)),
             annot.by = c("orig.ident"),
             cluster_cols = TRUE,
             fontsize = 7)

output <- getSignificance(ES2, group = "orig.ident", fit = "linear.model")


#Did not use at the end

df2 <- df1 %>% dplyr::select(1:50)
ES3 <- as.matrix(ES2)


Heatmap(t(scale(df2)), 
        name = "ssGSEA", #title of legend
        column_title = "ssGSEA",
        show_row_names = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 7), 
        col = colorRamp2(c(-3,0,3), c('blue', "white", "red")
        ))


ES2 <- ES %>% dplyr::select(1:20)
ES3 <- as.matrix(ES2)
#heatmap(ES3)


png("ESheatmap.png", height = 25 , width = 27, units="in", res=150)
Heatmap(scale(ES3), 
        name = "ssGSEA", #title of legend
        column_title = "ssGSEA",
        show_row_names = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 7), 
        col = colorRamp2(c(-3,0,3), c('blue', "white", "red"))
)

Heatmap(t(ES3), 
        name = "ssGSEA", #title of legend
        column_title = "ssGSEA",
        show_row_names = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 7), 
        col = colorRamp2(c(-3,0,3), c('blue', "white", "red"))
)

mel1 <- ES[grep("Melanoma", rownames(ES)), ] #18474x4872
skin <- ES[grep("Skin", rownames(ES)), ] #7122x4872
nMcSC <- ES[grep("nMcSC", rownames(ES)), ] #1430x4872
Tbpt <- ES[grep("Tbpt", rownames(ES)), ] #2965x4872

mel1 <- as.data.frame(mel1)
skin <- as.data.frame(skin)
nMcSC <- as.data.frame(nMcSC)
Tbpt <- as.data.frame(Tbpt)

df1 <- data.frame(rbind(df,df,df,df))
Mel2 <-colMeans(mel1)
Mel2 <- as.data.frame.list(Mel2)
skin2 <-colMeans(skin)
skin2 <- as.data.frame.list(skin2)
nMcSC2 <-colMeans(nMcSC)
nMcSC2 <- as.data.frame.list(nMcSC2)
Tbpt2 <-colMeans(Tbpt)
Tbpt2 <- as.data.frame.list(Tbpt2)
df1 <- data.frame(rbind(Mel2,skin2,nMcSC2,Tbpt2))
rownames(df1) <- c( "Melanoma", "skin", "nMcSC", "Tbpt")