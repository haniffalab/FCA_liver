args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))
if(length(args) != 8){
  stop('This pipeline requires 8 parameters: organ, ProjectName, save.at (name of RDS file where processed data are saved), sequencing.types (normal, 5GEX or VDJ), annotate.cells (boolean), identify.doublets (boolean), cell.type.SVM (folder name where cell type svm classifier for given organ is stored), doublet.svm (folder name where singlet/doublet svm classifier for given organ is stored);')
}

arguments.list = "
organ.arg             = args[1]
ProjectName.arg       = args[2]
save.at.arg           = args[3]
sequencing.types.arg  = args[4]
annotate.cells.arg    = args[5]
identify.doublets.arg = args[6]
cell.type.SVM.arg     = args[7]
doublet.svm.arg       = args[8]
"

eval(parse(text = arguments.list))

arguments.list = unlist(strsplit(arguments.list, "\n"))
arguments.list = arguments.list[!(arguments.list == "")]

for(n in 1:length(arguments.list)){
  argument = arguments.list[n]
  argument = gsub(pattern=" ", replacement="", x=argument)
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  argument.content = eval(parse(text = argument.name))
  eval(parse(text = argument.content))
  if (!exists(variable.name)){
    stop(sprintf("Argument %s not passed. Stopping ... ", variable.name))
  }
}

# create required folders for output and work material
output_folder = gsub(pattern="^\\d+_", replacement="", x=basename(getwd()))
output_folder = paste(output_folder, save.at, sep = "_")
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

source("../../tools/bunddle_utils.R")

sequencing.types = unlist(strsplit(sequencing.types, "-"))
save.at = file.path("../../data", save.at)
key.fname = "../../resources/key.csv"
sample.key.fname = "../../resources/sample_key.csv"

library(Seurat)
library(plyr)
library(dplyr)

# function to load sequentiall 10x data from many folders passed as a vector
# return a seurat object created from merging all the 10X data in the folders
# does not apply any filtering
load_10x_data_from_folders = function(folders, inside = "filtered/GRCh38/", key, sample.key){
  # load data from first folder. if only one folder then return
  sample.col.label = folders[1]
  nfolders = length(folders)
  prelabel.a = key$Prelabel[key$SUPPLIER.SAMPLE.NAME == sample.col.label]
  folder = key$SANGER.SAMPLE.ID[key$SUPPLIER.SAMPLE.NAME == sample.col.label]
  folder = file.path(folder, inside)
  print(sprintf("Loading sample %s from %s", prelabel.a, folder))
  seurat.obj.a = tryCatch({
    Read10X(folder)
  }, error = function(e){
    Read10X(file.path(unlist(strsplit(folder, '/'))[1], unlist(strsplit(folder, '/'))[2]))
  })
  colnames(seurat.obj.a) = paste(prelabel.a, colnames(seurat.obj.a), sep = "_")
  seurat.obj = CreateSeuratObject(raw.data = seurat.obj.a, min.cells = 0, min.genes = 0, project = "")
  # if there is only one folder to read data from than add the prelabel to the cell names and return the objec
  if (nfolders == 1){
    return(seurat.obj)
  }
  # if more the 1 folder load next folder(s)
  for (i in 2:nfolders){
    sample.col.label = folders[i]
    print(sample.col.label)
    prelabel.b = key$Prelabel[key$SUPPLIER.SAMPLE.NAME == sample.col.label]
    folder = key$SANGER.SAMPLE.ID[key$SUPPLIER.SAMPLE.NAME == sample.col.label]
    folder = file.path(folder, inside)
    print(sprintf("Loading sample %s from %s", prelabel.b, folder))
    seurat.obj.b = tryCatch({
       Read10X(folder)
    }, error = function(e){
      Read10X(file.path(unlist(strsplit(folder, '/'))[1], unlist(strsplit(folder, '/'))[2]))
    })
    seurat.obj = AddSamples(object=seurat.obj, new.data=seurat.obj.b, project=ProjectName, min.cells=0, min.genes=0, do.normalize=F, do.scale=F, do.center=F, add.cell.id=prelabel.b)
    print(seurat.obj)
  }
  # eliminate the multiple underscore heading from the cell names introduces by sequential merging
  cell.names = strsplit(colnames(seurat.obj@raw.data), "_")
  cell.names = lapply(cell.names, function(x)x[x != ""])
  cell.names = unlist(lapply(cell.names, FUN=function(parts){return(paste(parts, collapse="_"))}))
  colnames(seurat.obj@data) = cell.names
  colnames(seurat.obj@raw.data) = cell.names
  names(seurat.obj@ident) = cell.names
  return(seurat.obj)
}

# function to add meta data to a seurat object based on parsing cell names
# currently it adds: fetal ids, sort ids, tissue, lane, stage and sample type
add.meta.data = function(seurat.obj, sample.key, key){
  cell.names = strsplit(colnames(seurat.obj@data), "_")
  fetal.ids  = as.factor(unlist(lapply(cell.names, "[", 1)))
  tissue     = as.factor(unlist(lapply(cell.names, "[", 2)))
  sort.ids   = as.factor(unlist(lapply(cell.names, "[", 3)))
  lanes      = as.factor(unlist(lapply(cell.names, "[", 4)))
  key.key = which(sample.key$Sample %in% levels(fetal.ids))
  # map stages
  stages = plyr::mapvalues(x=fetal.ids, from=sample.key$Sample[key.key], to = sample.key$Stage[key.key])
  # map sample type
  sample.type = plyr::mapvalues(x=fetal.ids, from=sample.key$Sample[key.key], to = sample.key$Type[key.key])
  # map fetal ids
  fetal.ids = plyr::mapvalues(x=fetal.ids, from=sample.key$Sample[key.key], to = sample.key$Name[key.key])
  # create gender vector
  gender = strsplit(as.vector(fetal.ids), "_")
  gender = as.factor(unlist(lapply(gender, "[",2)))
  # create the AnnatomicalPart vector
  unique.lanes = as.vector(unique(lanes))
  unique.key = key[key$SANGER.SAMPLE.ID %in% unique.lanes, ]
  AnnatomicalPart = plyr::mapvalues(x=lanes, from=unique.key$SANGER.SAMPLE.ID, to=unique.key$AnnatomicalPart)
  # add the meta data
  seurat.obj@meta.data$fetal.ids       = fetal.ids
  seurat.obj@meta.data$sort.ids        = sort.ids 
  seurat.obj@meta.data$tissue          = tissue
  seurat.obj@meta.data$lanes           = lanes
  seurat.obj@meta.data$stages          = stages
  seurat.obj@meta.data$sample.type     = sample.type 
  seurat.obj@meta.data$gender          = gender
  seurat.obj@meta.data$AnnatomicalPart = AnnatomicalPart 
  return(seurat.obj)
}

# function to perform filtering on a seurat object
# this ensures all the datasets in a project are filtered with same criteria
filter.seurat = function(seurat.obj, min.cells = 3, min.genes = 200, project.name = "", mito.genes.treshold = .2){
  # apply filtering based on min.genes and min.cells
  print("Filtering on cell and gene numbers ... ")
  seurat.obj = CreateSeuratObject(raw.data = seurat.obj@raw.data, min.cells = min.cells, 
                                   min.genes = min.genes, project = "")
  seurat.obj.meta.data = seurat.obj@meta.data
  saveRDS(seurat.obj.meta.data, file.path(output_folder, 'meta_data_mingenes.RDS'))
  # calculate percentage of mitocondrial genes
  mito.genes = grep(pattern = "^MT-", x = rownames(x = seurat.obj@data), value = TRUE)
  percent.mito = Matrix::colSums(seurat.obj@raw.data[mito.genes, ])/Matrix::colSums(seurat.obj@raw.data)
  seurat.obj = AddMetaData(object = seurat.obj, metadata = percent.mito, col.name = "percent.mito")
  # filter on mitocondrial genes > mito.genes.treshold
  print("Filtering on mitochondrial genes")
  seurat.obj = FilterCells(object = seurat.obj, subset.names = c("percent.mito"), low.thresholds = c(-Inf), 
                            high.thresholds = c(mito.genes.treshold))
  seurat.obj.meta.data = seurat.obj@meta.data
  saveRDS(seurat.obj.meta.data, file.path(output_folder, 'meta_data_mitogenes.RDS'))
  return(seurat.obj)
}

# load the key
# then do View(key) to look for the datasets you required
# write the names of interest from key$V2 into data.folders
# you can now access the folder names that need to be uploaded for creating a seurat object with required data
#key = read.csv(file = key.fname, stringsAsFactors = FALSE, header=T)
key = read.csv(file = key.fname, stringsAsFactors = FALSE, header=T, sep="\t")
key$Fetus = unlist(regmatches(key$SUPPLIER.SAMPLE.NAME, gregexpr(pattern="F[0-9]{2}", text=key$SUPPLIER.SAMPLE.NAME)))
key = key[key$Organ != "other", ]


# load sample key
sample.key = read.csv(sample.key.fname, stringsAsFactors = F, sep = "\t")
# check all sample names in the key are also in the sample key:
print(paste("All sample names in the key are also in the sample_key: ", all(key$Fetus %in% sample.key$Sample), sep = ""))

# create prelabel column in key data frame
# the prelabel is attached to the cell names before each barcode
prelabel = paste(key$Fetus, key$Organ, sep = "_")
# the first 20 supplier labels have 4 fields so the gate field is at position 4
gate = strsplit(key$SUPPLIER.SAMPLE.NAME, split="_")
gate = as.vector(unlist(lapply(gate, "[", 3)))
gate = plyr::mapvalues(x = gate, from = c("CD45P", "CD45N", "TOT"), 
                        to = c("CD45+", "CD45-", "Total"))
prelabel = paste(prelabel, gate, sep = "_")
prelabel = paste(prelabel, key$SANGER.SAMPLE.ID, sep = "_")
key$Prelabel = prelabel

# next the data can be parsed by tissue

##########################################################################################
##########################################################################################
##########################################################################################

data.folders = key$SUPPLIER.SAMPLE.NAME[(key$Organ == organ & key$Sequencing %in% sequencing.types ) & key$Passed]
key = key[(key$Organ == organ & key$Sequencing %in% sequencing.types ) & key$Passed, ]
print("Next is the key: ")
print(key)

# load the data
cur.dir = getwd()
setwd("../../data/sc_count_matrices/")
seurat.obj = load_10x_data_from_folders(folders=data.folders, key = key, sample.key = sample.key)
setwd(cur.dir)

# parse meta data from cell names
seurat.obj = add.meta.data(seurat.obj, sample.key = sample.key, key = key)
print("Number of cells by lanes and gates before filtering:")
print(table(seurat.obj@meta.data$fetal.ids, seurat.obj@meta.data$sort.ids))

print('Cells by samples and gates before filtering:')
print(table(seurat.obj@meta.data$fetal.ids, seurat.obj@meta.data$sort.ids))

# apply filtering
seurat.obj = filter.seurat(seurat.obj=seurat.obj, project.name="")
# parse meta data from cell names because meta.data is lost during filtering
seurat.obj = add.meta.data(seurat.obj, sample.key = sample.key, key=key)

print('Cells by samples and gates after filtering:')
print(table(seurat.obj@meta.data$fetal.ids, seurat.obj@meta.data$sort.ids))

print('Cells by lanes and gates after filtering:')
print(table(seurat.obj@meta.data$lanes, seurat.obj@meta.data$sort.ids))

# normaliza data
print("Normalizing data ...")
seurat.obj = NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

print("Computing variable genes ...")
# find variable genes
seurat.obj = FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                                        dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                        x.high.cutoff = 3, y.cutoff = .625)
# calculate percentage of variable genes
print(paste("Percentage of variable genes:", round(100 * length(seurat.obj@var.genes) / dim(seurat.obj@data)[1], digits = 2), sep = " "))

# scale data in variable genes, otherwise pca is not possible
print("Scaling data ...")
seurat.obj = ScaleData(object=seurat.obj)

# run PCA
print("Performing PCA ...")
seurat.obj = RunPCA(object = seurat.obj, pc.genes = seurat.obj@var.genes, do.print = FALSE, pcs.print = 1:20, genes.print = 10)

# run TSNE
print("Performing TSNE")
seurat.obj = RunTSNE(object=seurat.obj, dims.use=1:20, seed.use=42, do.fast=TRUE)

# run umap
print("running UMAP")
umap.coordinates = RunUMAP(pca.df=seurat.obj@dr$pca@cell.embeddings, tool_addr=tool_addr, python.addr=python.addr)
rownames(umap.coordinates) = names(seurat.obj@ident)
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="cell.embeddings", new.data=as.matrix(umap.coordinates))
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="key", new.data="umap")

# run force-directed graph
print("Running force directed graph")
seurat.obj = BuildSNN(object=seurat.obj, reduction.type="pca", dims.use=1:20, plot.SNN=F)
fdg_coordinates = runFDG(pca.df=seurat.obj@dr$pca@cell.embeddings, snn=seurat.obj@snn, iterations=600, tool_addr=tool_addr, python.addr=python.addr)
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot="cell.embeddings", new.data=as.matrix(fdg_coordinates))
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot = "key", new.data = "fdg")

if(annotate.cells){
  print("Annotating cells ... ")
  seurat.obj@meta.data$cell.labels = Apply_Classifier_On_Seurat_Object(seurat.obj=seurat.obj, classifier.fname=cell.type.SVM, tool_addr=tool_addr, python.addr=python.addr)
}

if (identify.doublets){
  print("identifying doublets")
  seurat.obj@meta.data$doublets = Apply_Classifier_On_Seurat_Object(seurat.obj=seurat.obj, classifier.fname=doublet.svm, tool_addr=tool_addr, python.addr=python.addr)
}

print("saving data")
saveRDS(seurat.obj, save.at)


if (identify.doublets){
  print("Doublets and singlets: ")
  print(table(seurat.obj@meta.data$fetal.ids, seurat.obj@meta.data$doublets)) 
}

print("Ended beautifully ... ")

