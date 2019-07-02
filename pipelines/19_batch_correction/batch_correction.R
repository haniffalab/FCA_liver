args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg = args[1]
correct.by.arg  = args[2]
save.at.arg     = args[3]
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
save.at     = file.path("../../data", save.at)

source("../../tools/bunddle_utils.R")

library(harmony)
library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(methods)
library(utils)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

print("Performing PCA")
seurat.obj = RunPCA(object = seurat.obj, pc.genes = seurat.obj@var.genes, do.print = FALSE)

# run harmany batch correction
seurat.obj = RunHarmony(object=seurat.obj, group.by.vars=correct.by, theta=3)

print("Performing tSNE")
seurat.obj = RunTSNE(object=seurat.obj, dims.use=1:20, seed.use=42, do.fast=TRUE, reduction.use="harmony")

# compute UMAP based on harmomy results
print("running UMAP")
umap.coordinates = RunUMAP(pca.df=seurat.obj@dr$harmony@cell.embeddings, tool_addr=tool_addr, python.addr=python.addr)
rownames(umap.coordinates) = names(seurat.obj@ident)
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="cell.embeddings", new.data=as.matrix(umap.coordinates))
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="key", new.data="umap")

# run force-directed graph
print("Running force directed graph")
seurat.obj = BuildSNN(object=seurat.obj, reduction.type="harmony", dims.use=1:20, plot.SNN=F, force.recalc=T)
fdg_coordinates = runFDG(pca.df=seurat.obj@dr$harmony@cell.embeddings, snn=seurat.obj@snn, iterations=2000, tool_addr=tool_addr, python.addr=python.addr)
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot="cell.embeddings", new.data=as.matrix(fdg_coordinates))
seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot = "key", new.data = "fdg")

# saving the result
print("Saving data ...")
saveRDS(seurat.obj, save.at)

print("Ended beautifully ... ")
