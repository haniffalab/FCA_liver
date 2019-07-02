args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg  = args[1]
do.normalize.arg = args[2]
add.PCA.arg      = args[3]
add.TSNE.arg     = args[4]
add.UMAP.arg     = args[5]
add.FDG.arg      = args[6]
save.dr.arg      = args[7]
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

# create required folders for output and work material
output_folder = gsub(pattern="^\\d+_", replacement="", x=basename(getwd()))
output_folder = paste(output_folder, seurat.addr, sep = "_")
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

seurat.addr = file.path("../../data", seurat.addr)

source("../../tools/bunddle_utils.R")

library(Seurat)
library("sva")
library(plyr)
library(dplyr)
library(reshape)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

if(do.normalize){
  print("Normalizing data ... ")
  seurat.obj = NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  print("Applying COMBAT ...")
  expression.data = as.matrix(seurat.obj@data[seurat.obj@var.genes, ])
  pheno.data = data.frame(sample = names(seurat.obj@ident), 
                          batch = seurat.obj@meta.data$fetal.ids, 
                          stages = seurat.obj@meta.data$stages)
  batch = as.numeric(pheno.data$batch)
  pheno.data$batch = batch
  mod = model.matrix(~as.factor(stages), data=pheno.data)
  colnames(expression.data) = gsub(pattern="CD45[+]", replacement="CD45Pos", x=colnames(expression.data))
  colnames(expression.data) = gsub(pattern="CD45[-]", replacement="CD45Neg", x=colnames(expression.data))
  
  write.csv(expression.data, file.path(output_folder, "data.csv"), row.names = T)
  batch = data.frame(Batch = batch)
  rownames(batch) = colnames(expression.data)
  write.csv(batch, file.path(output_folder, "batch.csv"))
  
  command = sprintf("%s combat.py %s", python.addr, output_folder)
  system(command, wait = T)
  rm(mod, expression.data)
  
  combat_data = read.csv(file.path(output_folder, "combat.csv"), sep = ",", row.names = 1)
  print("COMBAT data loaded.")
  
  colnames(combat_data) = gsub(pattern="CD45Pos", replacement="CD45+", x=colnames(combat_data))
  colnames(combat_data) = gsub(pattern="CD45Neg", replacement="CD45-", x = colnames(combat_data))

  combat_data = as(as.matrix(combat_data), "dgCMatrix")
  genes.not = rownames(seurat.obj@data)[!(rownames(seurat.obj@data) %in% rownames(combat_data))]
  all.expression = seurat.obj@data[genes.not, ]
  all.expression = rbind(all.expression, combat_data)
  all.expression = all.expression[rownames(seurat.obj@data), ]
  seurat.obj@data = combat_data
  
  print("Scaling data ... ")
  seurat.obj = ScaleData(object=seurat.obj)
}

if(add.PCA){
  print("Performing PCA")
  seurat.obj = RunPCA(object = seurat.obj, pc.genes = seurat.obj@var.genes, do.print = FALSE)
}

if (add.TSNE){
  print("Performing tSNE")
  seurat.obj = RunTSNE(object=seurat.obj, dims.use=1:20, seed.use=42, do.fast=TRUE)
}

if (add.UMAP){
  # run umap
  print("running UMAP")
  umap.coordinates = RunUMAP(pca.df=seurat.obj@dr$pca@cell.embeddings, tool_addr=tool_addr, python.addr=python.addr)
  rownames(umap.coordinates) = names(seurat.obj@ident)
  seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="cell.embeddings", new.data=as.matrix(umap.coordinates))
  seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="key", new.data="umap")
}

if (add.FDG){
  # run force-directed graph
  print("Running force directed graph")
  seurat.obj = BuildSNN(object=seurat.obj, reduction.type="pca", dims.use=1:20, plot.SNN=F, force.recalc=T)
  fdg_coordinates = runFDG(pca.df=seurat.obj@dr$pca@cell.embeddings, snn=seurat.obj@snn, iterations=2000, tool_addr=tool_addr, python.addr=python.addr)
  seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot="cell.embeddings", new.data=as.matrix(fdg_coordinates))
  seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot = "key", new.data = "fdg")
}

print("Saving Seurat object")
saveRDS(seurat.obj, seurat.addr)

if(save.dr){
  CellNames  = as.vector(names(seurat.obj@ident))
  tSNEdata   = seurat.obj@dr$tsne@cell.embeddings
  UMAPdata   = seurat.obj@dr$umap@cell.embeddings
  FDGdata    = seurat.obj@dr$fdg@cell.embeddings
  PCAdata    = seurat.obj@dr$pca@cell.embeddings
  colnames(tSNEdata) = c("tSNEx", "tSNEy")
  colnames(UMAPdata) = c("UMAPx", "UMAPy")
  colnames(FDGdata)  = c("FDGx", "FDGy")
  dr_md_df = data.frame(CellNames = CellNames)
  dr_md_df = cbind(dr_md_df, tSNEdata, UMAPdata, FDGdata, PCAdata, seurat.obj@meta.data)
  save.to = file.path(output_folder, "dr_and_metadata.csv")
  write.csv(dr_md_df, save.to)
}else{
  unlink(output_folder, recursive=T, force=T) 
}

file.remove("Rplots.pdf")

print("Ended beautifully ... ")
