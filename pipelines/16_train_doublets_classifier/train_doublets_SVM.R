args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg = args[1]
save.to.arg     = args[2]
"

expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lapply(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
  stop(sprintf('This pipeline requires %s parameters: '))
}

eval(parse(text = arguments.list))

for(n in 1:length(expected_arguments)){
  argument = expected_arguments[n]
  argument = gsub(pattern=" ", replacement="", x=argument)
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  argument.content = eval(parse(text = argument.name))
  eval(parse(text = argument.content))
  if (!exists(variable.name)){
    stop(sprintf("Argument %s not passed. Stopping ... ", variable.name))
  }
}

seurat.addr = file.path("../../data", seurat.addr)

source("../../tools/bunddle_utils.R")

library(Seurat)
library(RColorBrewer)
library(plyr)
library(dplyr)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

save.to = file.path("../../resources", save.to)

dir.create(save.to)

# for HdCA samples:
#seurat.obj@meta.data$lanes <- seurat.obj@meta.data$biomaterial_id

# mix cells within lanes
lanes <- as.vector(unique(seurat.obj@meta.data$lanes))
mixersA <- c(); mixersB <- c();
for (i in 1:length(lanes)){
  lane = lanes[i]
  print(lane)
  lane.indices <- which(seurat.obj@meta.data$lanes == lane)
  mixersA <- c(mixersA, sample(x=lane.indices, size=length(lane.indices), replace=T))
  mixersB <- c(mixersB, sample(x=lane.indices, size=length(lane.indices), replace=T))
}

# create a matrix of raw data for singlets - these are the cell collected
singlets <- seurat.obj@raw.data[rownames(seurat.obj@data), colnames(seurat.obj@data)]

# create doublets: 1) select from the singlets matrix 2 other matrices using the mixersA and mixersB; 2) add these 2 matrices
doublets <-  singlets[, mixersA] + singlets[, mixersB]

# make the cell names unique
colnames(singlets) <- paste("Singlet",  1:dim(singlets)[2], sep = "_")
colnames(doublets) <- paste("Doublets", 1:dim(doublets)[2], sep = "_")

# merge the singlets and doublets sparse matrices
# need to look into whether cBind deprecation will affect results
merged.matrix <- Matrix::cBind(singlets, doublets)

# free up some space
rm(seurat.obj, lanes, mixersA, mixersB)

# create label vector
labels <- rep(c("Singlet", "Doublet"), each = ncol(singlets))
labels <- data.frame(Labels = labels)

# create a seurat object from singlets and doublets
seurat.obj <- CreateSeuratObject(raw.data = merged.matrix, project="Doublet_Classifier", min.cells=0, min.genes=0)

print("Checkpoint 2")
# free up some space
rm(singlets, doublets, merged.matrix)
print("Checkpoint 3")
# normalize the data
seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# compute variable genes
seurat.obj <- FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                                dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                x.high.cutoff = 3, y.cutoff = .625)
# save variable genes in the same folder as the classifier
classifier.features <- seurat.obj@var.genes
saveRDS(classifier.features, file.path(save.to, "feature_genes.RDS"))
print("Checkpoint 4")

shuffle.key <- sample(x=1:dim(labels)[1], size=dim(labels)[1], replace=F)

# save the labels and normalized data to disk 
write.csv(labels[shuffle.key,], file.path(save.to, "labels.csv"), row.names = F)
x.data <- as.data.frame(t(as.matrix(seurat.obj@data[classifier.features, ])))[shuffle.key,]
print(dim(x.data))
write.csv(x.data, file.path(save.to, "data.csv"), row.names = T)

command = sprintf("%s svm.py %s", python.addr, save.to)
system(command, wait = T)

file.remove(file.path(save.to, "data.csv"), file.path(save.to, "labels.csv"), './Rplots.pdf')

print("Ended beautifully ... ")
