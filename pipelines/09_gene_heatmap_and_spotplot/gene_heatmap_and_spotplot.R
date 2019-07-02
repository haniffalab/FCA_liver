args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))
if(length(args) != 8){
  stop('This pipeline requires 8 parameters: seurat.addr\n set.ident \n genes.to.plot (name of file containing genes to plot)\n cell.types (name of file containing cell types to plot or all)\n cluster.genes (boolean)\n diagonalize (boolean indicating to compute the order of the genes in such a way as to make the appearance of a diagonal on the spotplot - will override the cluster.genes boolean)\n plot.dims (4 tuple for plots dimenssion)\n save.gene.order (boolean indicate to save the genes in the order computed by clustering and/or diaganolization - can be NA or a file name)')
}

arguments.list = "
seurat.addr.arg     = args[1]
set.ident.arg       = args[2]
genes.to.plot.arg   = args[3]
cell.types.arg      = args[4]  
cluster.genes.arg   = args[5]
diagonalize.arg     = args[6]
plot.dims.arg       = args[7]
save.gene.order.arg = args[8]
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
output_folder = paste("09_gene_expression_heatmap_and_spotplot", seurat.addr, sep = "_")
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

seurat.addr   = file.path("../../data", seurat.addr)
genes.to.plot = file.path("../../resources", genes.to.plot)

library(Seurat)
library(dplyr)
library(plyr)
library(reshape)

#######################################################################################################

# load seurat object
print("loading data ...")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

# set identies
seurat.obj = SetAllIdent(object=seurat.obj, id=set.ident)

# load genes
genes.to.plot = file(genes.to.plot, "r")
genes = readLines(genes.to.plot, warn=F)
close(genes.to.plot)
genes = unlist(strsplit(genes, "\n"))

# check that all genes are in the dataset
if(!(all(genes %in% rownames(seurat.obj@data)))){
  not.found = genes[!(genes %in% rownames(seurat.obj@data))]
  print(sprintf("The following genes were not found in the data: %s", paste(not.found, collapse = ", ")))
  genes = genes[genes %in% rownames(seurat.obj@data)]
}

# check for duplicates
if(length(genes) != length(unique(genes))){
  duplicates = names(table(genes)[table(genes) > 1])
  duplicates = paste(duplicates, collapse = ", ")
  print(sprintf("Duplicates found: %s", duplicates))
  print("This will not affect the workflow, but be aware the heat map will have fewer genes than expected.")
  genes = unique(genes)
}

# rearange expression matrix by the order in cell types
if(cell.types != "all"){
  cell.types = file.path("../../resources", cell.types)
  cell_types_file = file(cell.types, "r")
  cell.types      = readLines(cell_types_file, warn = F)
  close(cell_types_file)
  cell.types = unlist(strsplit(cell.types, ", "))
  print(cell.types)
  print("All cell types in data set:")
  print(table(cell.types %in% as.vector(unique(seurat.obj@ident))))
}else{
  cell.types = sort(as.vector(unique(seurat.obj@ident)))
}

# subset expression data matrix
keep.cell.names = names(seurat.obj@ident)[seurat.obj@ident %in% cell.types]
expression.data = data.matrix(seurat.obj@data[genes, keep.cell.names])

# create a data matrix with mean expression of each marker by cell type
expression.data = t(expression.data)
expression.data = as.data.frame(expression.data)
expression.data = cbind(data.frame(CellLabels = as.vector(seurat.obj@ident[keep.cell.names])), expression.data)
expression.data = aggregate(expression.data[2:dim(expression.data)[2]], list(expression.data$CellLabels), mean)
expression.data = cbind(data.frame(CellType = expression.data$Group.1), expression.data[, 2:dim(expression.data)[2]])
rownames(expression.data) = expression.data$CellType
expression.data = expression.data[, 2:ncol(expression.data)]

# cluster the genes and reorder the expression matrix
if (cluster.genes){
  expression.distance = dist(x=t(expression.data), method="euclidian")
  gene.order = hclust(d=expression.distance, method="ward.D")$order
  expression.data = expression.data[, gene.order]
}

if (diagonalize){
  computer.vector.weight.center = function(vecn){
    indices = 1:length(vecn)
    sum(indices * vecn) / sum(vecn)
  }
  centers = apply(X=expression.data, MARGIN=2, FUN=computer.vector.weight.center)
  centers = order(centers)
  print(length(centers))
  expression.data = expression.data[, centers]
}

# plot the heatmap
expression.melt = reshape::melt(data=as.matrix(expression.data))
colnames(expression.melt) = c("CellTypes", "Genes", "Values")
expression.melt$CellTypes = factor(as.vector(expression.melt$CellTypes), levels = cell.types)
heatmap.plot = ggplot(expression.melt, aes(factor(Genes, levels = colnames(expression.data)), factor(CellTypes, levels = rev(cell.types)))) + geom_tile(aes(fill = Values), color = "black")
heatmap.plot = heatmap.plot + scale_fill_gradient(low = "lightblue", high = "darkred")
heatmap.plot = heatmap.plot + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                                     axis.text.x = element_text(angle = 45, hjust = 1))
heatmap.plot = heatmap.plot + labs(fill='Expression') 
heatmap.fname = file.path(output_folder, "./heatmap.pdf")
pdf(heatmap.fname, width = plot.dims[1], height = plot.dims[2])
print(heatmap.plot)
dev.off()

# plot diag matrix as dot plot (spot plot)
expression.data.r = expression.data
expression.data.r = expression.data.r[rev(cell.types), rev(colnames(expression.data.r))]
expression.melt = reshape::melt(data=as.matrix(expression.data.r))
colnames(expression.melt) = c("CellTypes", "Genes", "Values")
expression.melt$X = rep(1:length(unique(expression.melt$Genes)), each=nrow(expression.data.r))
expression.melt$Y = rep(length(unique(expression.melt$CellTypes)):1, times=ncol(expression.data.r))
colnames(expression.melt) = c("CellTypes", "Genes", "Expression", "X", "Y" )

max.expression = floor(max(expression.melt$Expression)) + 1

spot.plot = ggplot(expression.melt, aes(x = Y, y = X)) +
  geom_point(aes(size = Expression, color = Expression)) + 
  scale_color_gradient(limits = c(0, max.expression), breaks = seq(0,max.expression, by = 1), low = "lightsteelblue1", high = "darkred") +
  guides(color = guide_legend(), size = guide_legend()) +
  scale_size_continuous(limits=c(0, max.expression), breaks=seq(0, max.expression, by=1)) + 
  scale_y_discrete(name ="",  limits=colnames(expression.data.r)) +  
  scale_x_discrete(name ="",  limits=rev(rownames(expression.data.r))) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

splotplot.fname = file.path(output_folder, "./spotplot.pdf")
pdf(splotplot.fname, width = plot.dims[3], height = plot.dims[4])
print(spot.plot)
dev.off()

if(!is.na(save.gene.order)){
  save.gene.order = file.path("../../resources", save.gene.order)
  save.gene.order = file(save.gene.order, "w")
  writeLines(colnames(expression.data), save.gene.order)
  close(save.gene.order)
}

print("Ended beautifully ... ")
