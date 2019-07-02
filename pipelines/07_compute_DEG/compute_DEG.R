# seurat object address
seurat.addrs  <- "../../data/dummydata.RDS"
save.to       <- "../../resources/dummydata_markergenes.csv"
DE.downsample <- F  # flag to indicate to downsample by cluster before computing DE genes
category      <- "Annotation_5" # categories to calculate DE genes for

library(Seurat)
library(dplyr)

# load the seurat object
print("Loading seurat object ...")
seurat.obj <- readRDS(seurat.addrs)
seurat.obj <- SetAllIdent(object=seurat.obj, id=category)

# writing marker genes to disk
if (DE.downsample){
  cluster.ids <-unique(as.vector(seurat.obj@ident))
  cells.to.keep <- c()
  for (k in 1:length(cluster.ids)){
    cluster.id <- cluster.ids[k]
    cell.ids <- names(seurat.obj@ident)[seurat.obj@ident == cluster.id]
    cell.ids <- which(names(seurat.obj@ident) %in% cell.ids )
    cells.to.keep <- c(sample(x=cell.ids, size=min(200, length(cell.ids)), replace=F), cells.to.keep)
  }
  seurat.obj_d <- SubsetData(object=seurat.obj, cells.use=names(seurat.obj@ident)[cells.to.keep])
  seurat.obj_d <- NormalizeData(object = seurat.obj_d, normalization.method = "LogNormalize", scale.factor = 10000)
  print("Calculating marker genes: finished subseting, currently actually calculating the markers ... ")
  marker.genes <- FindAllMarkers(object = seurat.obj_d, only.pos = F, min.pct = 0.25, thresh.use = 0.25,
                                 genes.use = rownames(seurat.obj@data), test.use = "wilcox",
                                 random.seed = 42, print.bar=T, do.print=T, max.cells.per.ident = 200)
}else{
  print("Calculating marker genes ...  ")
  marker.genes <- FindAllMarkers(object = seurat.obj, only.pos = F, min.pct = 0.25, thresh.use = 0.25,
                                 genes.use = rownames(seurat.obj@data), test.use = "wilcox",
                                 random.seed = 42, print.bar=T, do.print=T, max.cells.per.ident = 200)
}

print("Saving marker genes ... ")
print(save.to)
write.csv(marker.genes, save.to)

print("Finished!")
