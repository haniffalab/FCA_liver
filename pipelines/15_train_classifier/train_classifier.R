args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg  = args[1]
marker.genes.addr = args[2]
save.at          = args[3]
classifier       = args[4]
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
working_dir = paste(sample(LETTERS, 50, replace=T),collapse = '')
material_dir = file.path(working_dir, 'material')
output_dir = file.path(working_dir, 'output')
dir.create(working_dir)
dir.create(material_dir)
dir.create(output_dir)

save.at = file.path('../../resources', save.at)

seurat.addr = file.path("../../data", seurat.addr)

classifier = paste(classifier, '.py', sep = '.')

source("../../tools/bunddle_utils.R")

library(Seurat)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

# create and save label data frame
print("Create and save label data frame ...")
singlets <- as.vector(seurat.obj@meta.data$doublets) == "Singlet"
labels   <- data.frame(Labels = as.vector(seurat.obj@meta.data$cell.labels)[singlets])
write.csv(labels, file.path(material_dir, 'labels.csv'), row.names = F)

# save variable genes in the output folder
print("Choose features genes ...")
marker.genes = file.path('../../resources/marker_genes', marker.genes)
marker.genes <- read.csv(marker.genes)
marker.genes <- marker.genes %>% group_by(cluster) %>% top_n(20, avg_logFC)
classifier.features <- unique(as.vector(marker.genes$gene))
saveRDS(classifier.features, file.path(output_dir, 'feature_genes.RDS'))

# save the normalized data to disk
print("saving training data to disk ...")
cell.names <- names(seurat.obj@ident)[singlets]
x.data <- as.data.frame(t(as.matrix(seurat.obj@data[classifier.features, cell.names])))
write.csv(x.data, file.path(material_dir, 'data.csv'), row.names = T)

print("initiating SVM trainer ... ")
system(sprintf('%s svm.py %s %s', python.addr, material_dir, output_dir), wait = T)

# plot confusion matrix
cnf_matrix = read.csv(file.path(output_dir, 'confusion_matrix.csv'))
cnf_matrix <- cnf_matrix[, -c(1)]
confusion <- expand.grid(Actual = colnames(cnf_matrix), Predicted = colnames(cnf_matrix))
cnf_matrix <- cnf_matrix / colSums(cnf_matrix)
confusion$freq <- rapply(cnf_matrix, c)
pdf(file.path(output_dir, 'confusion_matrix.pdf'), width = 14, height = 14)
ggplot(data = confusion, aes(x = Actual, y = Predicted)) + geom_tile(aes(fill = freq)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

unlink(save.at, recursive=T, force=T)
dir.create(save.at)

file.rename(from=file.path(output_dir, 'feature_genes.RDS'), to=file.path(save.at, "feature_genes.RDS"))
file.rename(from=file.path(output_dir, 'classification_report.txt'), to=file.path(save.at, "classification_report.txt"))
file.rename(from=file.path(output_dir, 'confusion_matrix.csv'), to=file.path(save.at, "confusion_matrix.csv"))
file.rename(from=file.path(output_dir, 'confusion_matrix.pdf'), to=file.path(save.at, "confusion_matrix.pdf"))
file.rename(from=file.path(output_dir, 'model.pickle'), to=file.path(save.at, "model.pickle"))
file.rename(from=file.path(output_dir, 'pca.pickle'), to=file.path(save.at, "pca.pickle"))

unlink(working_dir, recursive=T, force=T)

print("Ended beautifully ... ")
