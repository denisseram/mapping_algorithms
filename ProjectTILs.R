install.packages("remotes")
library(remotes)

remotes::install_github("carmonalab/STACAS")
remotes::install_github("carmonalab/ProjecTILs")

library(ProjecTILs)
library(Seurat)

ref <- load.reference.map()
data(query_example_seurat)

query.projected <- Run.ProjecTILs(query_example_seurat, ref=ref)

#Let's explore the reference dataset
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
DimPlot(ref, label = T, cols = refCols)

markers <- c("Cd4", "Cd8a", "Ccr7", "Tcf7", "Pdcd1", "Havcr2", "Tox", "Izumo1r",
             "Cxcr6", "Xcl1", "Gzmb", "Gzmk", "Ifng", "Foxp3")
VlnPlot(ref, features = markers, stack = T, flip = T, fill.by = "ident", cols = refCols,
        assay = "RNA") + NoLegend()

querydata <- ProjecTILs::query_example_seurat

library(GEOquery)
geo_acc <- "GSE86028"
getGEOSuppFiles(geo_acc)

fname2 <- sprintf("%s/GSE86028_TILs_sc_wt_mtko.tpm.log2.txt.gz", geo_acc)
querydata2 <- read.sc.query(fname2, type = "raw.log2")

query.projected <- Run.ProjecTILs(querydata, ref = ref)

plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5)


plot.statepred.composition(ref, query.projected, metric = "Percent")

genes4radar = c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Sell", "Gzmb", "Gzmk", "Pdcd1",
                "Havcr2", "Tox", "Mki67")

plot.states.radar(ref, query = query.projected, genes4radar = genes4radar, min.cells = 20)

