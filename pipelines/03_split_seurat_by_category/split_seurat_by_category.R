library(Seurat)

# can split RDS by doublets, or any other metadata column in the seurat object

sort.by <- "doublets"
seurat.addrs <- "/home/b8058304/single_cell_data_analysis_bundle/data/test_yolk_sac_all.RDS"

#######################################################################################################
#######################################################################################################
#######################################################################################################

# load the seurat object
print("Loading the data ... ")  
seurat.obj <- readRDS(seurat.addrs)

# get categories of sort.by
eval(parse(text=paste("cats <- seurat.obj@meta.data$", sort.by, sep = "")))
cats.unique <- unique(as.vector(cats))
print(cats.unique)

# set the ident slot to the sort.by meta.data
seurat.obj <- SetAllIdent(object=seurat.obj, id=sort.by)

# init list to store seurat objects splitted by cats
store <- c()
for (i in 1:length(cats.unique)){
  category <- cats.unique[i]
  cat.object <- SubsetData(object=seurat.obj, ident.use=category, subset.raw=T, do.clean = T)
  # filter out cells you do not want to keep
  cells.to.keep <- names(cat.object@ident)[cat.object@meta.data[sort.by] == category]
  cat.object <- SubsetData(object=cat.object, cells.use=cells.to.keep)
  
  print(paste("Doing: ", category, sep = ""))
  
  print("Check point 3")
  print(dim(cat.object@data))
  print(dim(cat.object@raw.data))
  
  file.name <- gsub(pattern=" ", replacement="_", x=category)
  file.name <- paste(file.name, ".RDS", sep = "")
  saveRDS(cat.object, file.name)
}

print("Ended beautifully")
