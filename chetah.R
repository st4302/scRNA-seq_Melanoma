###CHETAH###

library(scater)
library(Seurat)
library(cowplot)


#Convert to single cell experiment
#pbmc.sce <- as.SingleCellExperiment(dat)
#p2 <- plotPCA(pbmc.sce, colour_by = "ident")
#p2

Mouseseq <- MouseRNAseqData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))
celltypes_hn <- Mouseseq$label.main
counts_hn <- assay(Mouseseq)
class(counts_hn)
counts_hn[1:5,1:5]

counts <- GetAssayData(object=dat)
pbmc.sce <- as.SingleCellExperiment(dat)
umap <- reducedDim(pbmc.sce,type="UMAP")
input_mel <- SingleCellExperiment(assays = list(counts = counts),
                                  reducedDims = SimpleList(UMAP = umap))


class(counts)
umap[1:5,]
all.equal(rownames(umap), colnames(counts))
headneck_ref <- SingleCellExperiment(assays = list(counts = counts_hn),
                                     colData = DataFrame(celltypes = celltypes_hn))
input_mel <- CHETAHclassifier(input = input_mel,
                              ref_cells = headneck_ref)

PlotCHETAH(input = input_mel)
PlotCHETAH(input = input_mel, interm = TRUE)

input_mel <- Classify(input_mel, 0)
PlotCHETAH(input = input_mel, tree = FALSE)
