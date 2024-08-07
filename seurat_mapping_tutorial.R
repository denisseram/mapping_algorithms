library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)

# We use the data from seurat (Pancreas dataset)

InstallData("panc8")
panc8 <- LoadData("panc8")
#it comes from different technologies
table(panc8$tech)


# we will use data from 2 technologies for the reference
pancreas.ref <- subset(panc8, tech %in% c("celseq2", "smartseq2"))
pancreas.ref[["RNA"]] <- split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)

# pre-process dataset (without integration)
pancreas.ref <- NormalizeData(pancreas.ref)
pancreas.ref <- FindVariableFeatures(pancreas.ref)
pancreas.ref <- ScaleData(pancreas.ref)
pancreas.ref <- RunPCA(pancreas.ref)
pancreas.ref <- FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)

pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30)
DimPlot(pancreas.ref, group.by = c("celltytpe", "tech"))

#Now, they integrate the datasets into a shared reference
pancreas.ref <- IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.cca", verbose = FALSE)
pancreas.ref <- FindNeighbors(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)

pancreas.ref <- RunUMAP(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

#we use the TransferData() function to classify the query cells based on reference data (a vector of reference cell type labels).
# select two technologies for the query datasets
pancreas.query <- subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query <- NormalizeData(pancreas.query)
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)


#evaluation how well is the algorithm performing
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
table(pancreas.query$predicted.id)

VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

pancreas.query <- TransferData(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                               refdata = list(celltype = "celltype"))
pancreas.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                                      new.reduction.name = "ref.pca")
pancreas.query <- ProjectUMAP(query = pancreas.query, query.reduction = "ref.pca", reference = pancreas.ref,
                              reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
