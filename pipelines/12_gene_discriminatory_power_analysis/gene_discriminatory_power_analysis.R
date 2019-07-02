args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg     = args[1]
feature_genes.arg   = args[2]
cell.types.arg      = args[3]
save.to.dir.arg     = args[4]
ident.set.arg       = args[5]
type.to.colours.arg = args[6]
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

library(Seurat)
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(RColorBrewer)

source("../../tools/bunddle_utils.R")

seurat.addr        = file.path("../../data", seurat.addr)
cell.types         = file.path('../../resources', cell.types)
feature_genes      = file.path('../../resources', feature_genes)
type.to.colours    = file.path("../../resources", type.to.colours)

cell.types.file = file(cell.types, "r")
cell.types = readLines(cell.types.file)
close(cell.types.file)

feature_genes.file = file(feature_genes, "r")
feature_genes = readLines(feature_genes.file)
close(feature_genes.file)

################################################################################################################

# a plotting function for indexed legend
plot.indexed.legend = function(label.vector, color.vector, ncols = 2, left.limit = 3.4, symbol.size = 8, text.size = 10, padH = 1, padV = 1, padRight = 0){
  if (length(label.vector) != length(color.vector)){
    stop("number of labels is different from number colors\nAdvice: learn to count!")
  }
  if (length(ncol) > length(label.vector)){
    stop("You cannot have more columns than labels\nSolution: Learn to count")
  }
  indices.vector = 1:length(label.vector)
  label.no = length(label.vector)
  nrows = ceiling(label.no / ncols)
  legend.frame = data.frame(X = rep(0, label.no), Y = rep(0, label.no), CS = color.vector, Txt = label.vector)
  legend.frame$X = rep(1:ncols, each=nrows)[1:nrow(legend.frame)]
  legend.frame$Y = rep(nrows:1, times = ncols)[1:nrow(legend.frame)]
  Xrange = range(legend.frame$X)
  Yrange = range(legend.frame$Y)
  plot.obj = ggplot(data = legend.frame, aes(x = X, y = Y))
  plot.obj = plot.obj + geom_point(size = symbol.size, colour = color.vector)
  plot.obj = plot.obj + scale_x_continuous(limits = c(Xrange[1] - padRight, Xrange[2] + padH))
  plot.obj = plot.obj + scale_y_continuous(limits = c(Yrange[1] - padV, Yrange[2] + padV))
  plot.obj = plot.obj + theme_void()
  
  plot.obj = plot.obj + annotate("text", x=legend.frame$X, y = legend.frame$Y, label = indices.vector, size = text.size)
  plot.obj = plot.obj + annotate("text", x=legend.frame$X+.1, y = legend.frame$Y, label=legend.frame$Txt, hjust = 0, size = text.size)
  return(plot.obj)
}

# plotting function for dimensionaly-reduced data to label population by a round indexed label
dr.plot = function(point.labels, dr1, dr2, dr1.name, dr2.name, no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2, use.cols = NULL, use.labels = NULL, limits = NULL, annotate.plot = T){
  df.dr = data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2)
  if(is.null(use.labels)){
    p.labels = sort(unique(as.vector(point.labels)))
  }
  else{
    p.labels = use.labels
  }
  df.dr$Cell.Labels = factor(df.dr$Cell.Labels, levels=p.labels)
  p.labels.medians = aggregate(df.dr[, 2:3], list(df.dr$Cell.Labels), median)
  df.dr$Cell.Labels = mapvalues(x = df.dr$Cell.Labels, from = p.labels, to = paste(1:length(p.labels), p.labels, sep = " "))
  if(is.null(use.cols)){
    set.seed(random_state)
    plt.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(p.labels)))
  }else{
    plt.colours = use.cols
  }
  index.map = 1:length(p.labels)
  plot.obj = ggplot(data = df.dr, aes(x = DR1, y = DR2, color = Cell.Labels))
  plot.obj = plot.obj + geom_point(size = pt.size)
  plot.obj = plot.obj + scale_color_manual(values=plt.colours)
  if(annotate.plot){
    plot.obj = plot.obj + geom_point(data=p.labels.medians,aes(x = DR1, y = DR2), colour = "gray", size = plt.lb.sz, fill = plt.colours, alpha = .5, pch = 21)
    plot.obj = plot.obj + annotate("text", x=p.labels.medians$DR1, y = p.labels.medians$DR2, label = index.map, size = txt.lb.size)
  }
  if (no.legend){
    plot.obj = plot.obj + theme(legend.position="none")
  }else{
    plot.obj = plot.obj + guides(color = guide_legend(override.aes = list(size=5)))
  }
  plot.obj = plot.obj + xlab(dr1.name) + ylab(dr2.name)
  if(!is.null(limits)){
    X0 = limits[1]; X1 = limits[2]; Y0 = limits[3]; Y1 = limits[4];
    plot.obj = plot.obj + scale_x_continuous(limits = c(X0, X1))
    plot.obj = plot.obj + scale_y_continuous(limits = c(Y0, Y1))
  }
  return(plot.obj)
}

################################################################################################################

# load seurat object
print("loading data ...")
seurat.obj = readRDS(seurat.addr)
print("Loaded data")

# set the clustering identity
seurat.obj = SetAllIdent(object=seurat.obj, id = ident.set)

# subset data on cell types
seurat.obj = SubsetData(object=seurat.obj, ident.use=cell.types)

# select on singlets
seurat.obj = SetAllIdent(object=seurat.obj, id = "doublets")
seurat.obj = SubsetData(object=seurat.obj, ident.use=c("Singlet"))
seurat.obj = SetAllIdent(object=seurat.obj, id = ident.set)

# normaliza data
print("Normalizing data ...")
seurat.obj = NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# check that all genes are in the dataset
print("check that genes are in the dataset")
if(!(all(feature_genes %in% rownames(seurat.obj@data)))){
  not.found = feature_genes[!(feature_genes %in% rownames(seurat.obj@data))]
  print(not.found)
}

# check for duplicates
print("check for duplicates")
if(length(feature_genes) != length(unique(feature_genes))){
  duplicates = names(table(feature_genes)[table(feature_genes) > 1])
  duplicates = paste(duplicates, collapse = ", ")
  print(sprintf("Duplicates found: %s", duplicates))
  print("This will not affect the workflow, but be aware the heat map will have a smaller genes than expected.")
  feature_genes = unique(feature_genes)
}

# create folder for saving the results
print("creating folders")
dir.create(save.to.dir)

# create folder to save working material
material_folder = file.path(save.to.dir, "material")
unlink(material_folder, recursive=T, force=T)
dir.create(material_folder)

# subsetting seurat object so that we do not get a 'problem too large' error 
seurat.obj = SetAllIdent(seurat.obj, id="cell.labels")
seurat.obj = SubsetData(seurat.obj, max.cells.per.ident = 1000)

# write the cluster labels to disk
if (!is.na(type.to.colours)){
  type.to.colour = read.csv(type.to.colours)
  filter.key   =  type.to.colour$CellTypes %in% as.vector(unique(seurat.obj@ident))
  cell.labels  = as.vector(type.to.colour$CellTypes[filter.key])
  cell.colours = as.vector(type.to.colour$Colours[filter.key])
}else{
  cell.labels  = sort(as.vector(unique(seurat.obj@ident)))
  cell.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
}

labels   = data.frame(Labels = as.vector(seurat.obj@ident))
labels$Colours = mapvalues(x=labels$Labels, from=cell.labels, to=cell.colours)
write.csv(labels, file.path(material_folder, "labels.csv"), row.names = F)

# write the feature data to disk
print("writing data.csv")
matrix <- as.matrix(seurat.obj@data)
feature_matrix <- subset(matrix, rownames(matrix) %in% feature_genes)
x.data <- as.data.frame(t(matrix))
write.csv(x.data, file.path(material_folder, "data.csv"), row.names = T)

# run the random forest classifier and get the confusion matrix
print("running random forest classifier")
command = paste(python.addr, sprintf("random_forest_classifier.py %s", save.to.dir), sep = " ")
system(command, wait = T)

print("plot the confusion matrix")

cnf_matrix = read.csv(file.path(material_folder, "confusion_matrix.csv"), check.names = F)
cnf_matrix = cnf_matrix[, -c(1)]
confusion = expand.grid(Actual = colnames(cnf_matrix), Predicted = colnames(cnf_matrix))
cnf_matrix_colSums = colSums(cnf_matrix)
cnf_matrix_colSums[cnf_matrix_colSums == 0] = 1.0
cnf_matrix_colSums_matrix = matrix(ncol = length(cnf_matrix_colSums), nrow = length(cnf_matrix_colSums))
cnf_matrix_colSums_matrix[] = cnf_matrix_colSums
cnf_matrix = cnf_matrix / t(cnf_matrix_colSums_matrix)
confusion$Frequency = rapply(cnf_matrix, c)
confusion$Actual = factor(as.vector(confusion$Actual), levels = cell.labels)
confusion$Predicted = factor(as.vector(confusion$Predicted), levels = rev(cell.labels))
confusion.plot = ggplot(data = confusion, aes(x = Actual, y = Predicted)) + geom_tile(aes(fill = Frequency)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_gradient(low = "lightblue", high = "darkred")
pdf(file.path(save.to.dir, "confusion_matrix.pdf"), width = 14, height = 14)
print(confusion.plot)
dev.off()

print("Ended beautifully.")
