args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg = args[1]
save.at.arg     = args[2]
doublet.svm.arg = args[3]
"
{
expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lapply(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
  stop(sprintf('This pipeline requires these parameters: seurat.addr ; save.at (name of RDS file where processed data are saved), doublet.svm (folder name where singlet/doublet svm classifier for given organ is stored)'))
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

output_folder_material = file.path(output_folder, "material")
dir.create(output_folder_material)

seurat.addr = file.path("../../data", seurat.addr)

source("../../tools/bunddle_utils.R")

save.at = file.path("../../data", save.at)

library(Seurat)
library(RColorBrewer)
library(dplyr)
library(plyr)
}
###############################################################################
# Load data
print("Loading data ...")
seurat.obj = readRDS(seurat.addr) # seurat.addr.arg
print("Data loaded.")

print("Identifying doublets")
seurat.obj@meta.data$doublets <- Apply_Classifier_On_Seurat_Object(
  seurat.obj = seurat.obj,
  classifier.fname = doublet.svm, # doublet.svm.arg
  tool_addr = tool_addr,
  python.addr = python.addr)

print("Doublets and singlets: ")
print(table(seurat.obj@meta.data$lanes, seurat.obj@meta.data$doublets))

print("Saving data")
saveRDS(seurat.obj, save.at) # save.at.arg

unlink(output_folder_material, recursive = T, force = T)
