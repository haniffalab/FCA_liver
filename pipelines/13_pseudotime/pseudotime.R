args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))
args = gsub(pattern = '@@', replacement = ' ', x = args)

arguments.list = "
seurat.addr.arg     = args[1]
set.ident.arg       = args[2]
cell.types.arg      = args[3]
root_cell_type.arg  = args[4]
var.genes.arg       = args[5]
type.to.colours.arg = args[6]
"

expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lapply(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
  stop(sprintf('This pipeline requires %s parameters: ', length(expected_arguments)))
}

eval(parse(text = arguments.list))

for(n in 1:length(expected_arguments)){
  argument = expected_arguments[n]
  #argument = gsub(pattern=" ", replacement="", x=argument)
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  variable.name = gsub(pattern=" ", replacement="", argument.name)
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

seurat.addr   = file.path("../../data", seurat.addr)

source("../../tools/bunddle_utils.R")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(monocle)
library(dplyr)
library(reshape2)

#######################################################################################################
###########

print("printing cell.types")
print(cell.types)
print("printing root_cell_type")
print(root_cell_type)

ma = function(arr, kernel = 50){
  res = arr
  n = 2 * kernel
  for(i in 1:length(arr)){
    start_index = max(1, i - kernel)
    stop_index = min(length(arr), i + kernel)
    res[i] = mean(arr[start_index:stop_index])
  }
  res
}

adaptive.moving_average = function(arr, kernel = 10, minim_kernel = 10, range.factor = 5){
  res = arr
  n = 2 * kernel
  for(i in 1:length(arr)){
    start_index = max(1, i - kernel)
    stop_index = min(length(arr), i + kernel)
    local_sd = sd(arr[start_index:stop_index])
    local_kernel = minim_kernel + round(range.factor / (local_sd + .1))
    start_index = max(1, i - local_kernel)
    stop_index = min(length(arr), i + local_kernel)
    res[i] = mean(arr[start_index:stop_index])
  }
  res
}

###########
#######################################################################################################

print("Loading data ...")
seurat.obj = readRDS(seurat.addr)
seurat.obj = SetAllIdent(object=seurat.obj, id=set.ident)
print("Data loaded.")

print("Subseting data on singlets and required cell populations")
if(cell.types == "all"){
  cell.types = as.vector(unique(seurat.obj@ident))
}
print(table(seurat.obj@ident))

print("Subseting data ...")
to.keep = names(seurat.obj@ident)[as.vector(seurat.obj@ident) %in% cell.types]
seurat.obj = SubsetData(object=seurat.obj, cells.use=to.keep)
seurat.obj@ident = factor(seurat.obj@ident, levels = cell.types)
print(table(seurat.obj@ident))

print("Writing data to disk ...")
# save raw data to disk
raw_data = seurat.obj@raw.data
raw_data = raw_data[rownames(seurat.obj@data), colnames(seurat.obj@data)]

# decomment the next lines if there is a list of genes that you need to exclude
to_exclude    = readRDS('fca_cellcycle_genes.RDS')
genes_to_keep = rownames(raw_data)
genes_to_keep = genes_to_keep[!(genes_to_keep %in% to_exclude)]
raw_data = raw_data[genes_to_keep, colnames(seurat.obj@data)]

writeMM(raw_data, file.path(output_folder_material, "raw_data.mtx"))
# save gene names
gene_names = rownames(raw_data)
write.csv(data.frame(Genes = gene_names), file.path(output_folder_material, "genenames.csv"))
# save cell names
cell_names = colnames(raw_data)
write.csv(data.frame(Cells = cell_names), file.path(output_folder_material, "cellnames.csv"))
# write cell labels to disk
write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)

print("Computing pseudotime using pdt.scanpy.py...")
# compute pseudotime in python scanpy
command = sprintf("%s pdt_scanpy.py %s %s", python.addr, root_cell_type, output_folder)
system(command, wait=T)

print("finished running .py")

# get cell labels and colours
if (!is.na(type.to.colours)){
  type.to.colours = file.path("../../resources", type.to.colours)
  type.to.colour = read.csv(type.to.colours)
  print("printing type.to.colour after it is loaded in")
  print(type.to.colour)
  print("printing seurat obj idents which the typetocol arg will be compared against in next lines")
  print(as.vector(unique(seurat.obj@ident)))
  filter.key   = type.to.colour$CellTypes %in% as.vector(unique(seurat.obj@ident))
  cell.labels  = as.vector(type.to.colour$CellTypes[filter.key])
  cell.colours = as.vector(type.to.colour$Colours[filter.key])
}else{
  cell.labels  = sort(as.vector(unique(seurat.obj@ident)))
  cell.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
}
print("printing cell.labels")
print(cell.labels)
print("printing cell.colours")
print(cell.colours)

# load pseudotime
print('reading pseudotime values')
pseudotime = read.csv(file.path(output_folder_material, "pseudotime.csv"), row.names = 1, header = F)
print("Are the cells in the same order in both pseudotime and seurat object? ")
print(all(rownames(pseudotime) == names(seurat.obj@ident)))
pseudotime$CellTypes = seurat.obj@ident

colnames(pseudotime) = c("Pseudotime", "CellType")
pseudotime$Color = mapvalues(x=pseudotime$CellType, from=cell.labels, to=cell.colours)
pseudotime$Color = factor(as.vector(pseudotime$Color), levels = cell.colours)
pseudotime$CellType = factor(as.vector(pseudotime$CellType), levels = cell.labels)
colnames(pseudotime) = c("Pseudotime", "Cell Type", "Color")
# making sure that there are no inf values in pdt column
#pseudotime["Pseudotime"][pseudotime["Pseudotime"] == "Inf"] <- 1

plot.density = ggplot(data = pseudotime, aes(x = Pseudotime, color = `Cell Type`, fill = `Cell Type`)) + geom_density(alpha = .7)
plot.density = plot.density + scale_x_continuous(position = "top", limits = c(.0, 1.0), expand = c(0.0, .0))
plot.density = plot.density + scale_color_manual(values = cell.colours)
plot.density = plot.density + scale_fill_manual(values = cell.colours)
plot.density = plot.density + theme(axis.title.y = element_blank(), 
                                     axis.text.y  = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.line.y  = element_blank(),
                                     axis.title.x = element_text(size = 25),
                                     legend.position = c(0, 1),
                                     legend.justification = c(0, 1))

print("printing cell.colours which is used for scale_colour_manual for plot.density plot")
print(cell.colours)

# compute diff genes
print("Computing var genes by cell type...")
cds                   = newCellDataSet(cellData = as.matrix(raw_data), phenoData=NULL, featureData=NULL, expressionFamily = negbinomial.size())
print("printing cds made using newCellDataSet function")
print(cds)
pData(cds)$Cluster    = as.vector(seurat.obj@ident)
print("printing cds after adding cluster to pdata")
print(cds)
print("running estimatesizefactors for cds")
cds                   = estimateSizeFactors(cds)
pData(cds)$Pseudotime = pseudotime$Pseudotime

if (is.na(var.genes)){
  var.genes.total = c()
  print('Computing variable genes ... ')
  for (j in 1:length(cell.labels)){
    print(sprintf("Choice %s out of %s ... ", as.character(j), as.character(length(cell.labels))))
    choices = pseudotime$`Cell Type` == cell.labels[j]
    var.genes = differentialGeneTest(cds[, choices], fullModelFormulaStr = "~sm.ns(Pseudotime)")
    var.genes = cbind(var.genes, data.frame(gene_id = rownames(var.genes)))
    var.genes.ch = var.genes %>% arrange(qval)
    var.genes.ch = as.vector(var.genes.ch$gene_id[1:100])
    var.genes.total = union(var.genes.total, var.genes.ch)
  }
  
  print("Computing var genes globally...")
  var.genes = differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  var.genes = cbind(var.genes, data.frame(gene_id = rownames(var.genes)))
  var.genes.ch = var.genes %>% arrange(qval)
  var.genes.ch = as.vector(var.genes.ch$gene_id[1:100])
  var.genes.total = union(var.genes.total, var.genes.ch)
  MT_genes = var.genes.total[grep("^MT-", x=var.genes.total, ignore.case=T)]
  var.genes.total = setdiff(var.genes.total, MT_genes)
}else{
  var.genes.file = file.path('../../resources', var.genes)
  var.genes.file = file(var.genes.file)
  var.genes.total = readLines(var.genes.file)
  var.genes.total = as.vector(unique(var.genes.total))
  var.genes.total = var.genes.total[var.genes.total != '']
  close(var.genes.file)
}

# saving the genes to disk

print("Heavy computing finished. Next saving to output...")

print("calculating var_gene_expression")

# cluster genes based on their min-max normalized values
var_gene_expression = as.matrix(seurat.obj@data[var.genes.total, order(pseudotime$Pseudotime)])
var_gene_expression = t(apply(var_gene_expression, 1, adaptive.moving_average, kernel = 15, minim_kernel = 1, range.factor=15))

# min-max normalization
var_gene_min = apply(var_gene_expression, 1, min)
var_gene_expression = var_gene_expression - var_gene_min
var_gene_genes_max = apply(var_gene_expression, 1, max)
var_gene_expression = var_gene_expression / var_gene_genes_max

print("clustering genes by level of expression")

# actual clustering of genes
d_matrix = as.dist(1.0 - cor(t(as.matrix(var_gene_expression)), method="spearman"))
genes_clust = hclust(d=d_matrix, method="ward.D2")
genes.in.order = var.genes.total[genes_clust$order] 

# plot min-max normalized expression
###################################################################################################
raw_data_genes = as.matrix(seurat.obj@data[rev(genes.in.order), order(pseudotime$Pseudotime)])
raw_data_genes = t(apply(raw_data_genes, 1, adaptive.moving_average, kernel = 15, minim_kernel = 1, range.factor=15))

# min-max normalization
raw_data_genes_min = apply(raw_data_genes, 1, min)
raw_data_genes = raw_data_genes - raw_data_genes_min
raw_data_genes_max = apply(raw_data_genes, 1, max)
raw_data_genes = raw_data_genes / raw_data_genes_max

print("group genes by pdt")

# group by pdt
pdt = range(pseudotime$Pseudotime)
pdt = seq(pdt[1], pdt[2], length.out=100)
pdt_data = c()
for (k in 1:nrow(raw_data_genes)){
  for(j in 1:length(pdt)){
    local_pdt = pdt[j]
    pdt_index = abs(pseudotime$Pseudotime[order(pseudotime$Pseudotime)] - local_pdt)
    pdt_index = which(pdt_index == min(pdt_index))
    pdt_data = c(pdt_data, raw_data_genes[k, pdt_index])
  }
}
pdt_data = matrix(data=pdt_data, nrow=nrow(raw_data_genes), byrow=T)
print("printing nrow pdt_data")
nrow(pdt_data)
print("printing ncol pdt_data")
ncol(pdt_data)
rownames(pdt_data) = rownames(raw_data_genes)
colnames(pdt_data) = paste("PDT", 1:100, sep = "")
#colnames(pdt_data) = paste("PDT", 1:ncol(pdt_data), sep = "")

# smooth a bit the pdt_data matrx
pdt_data = t(apply(pdt_data, 1, ma, kernel = 7))
pdt_data = pdt_data - apply(pdt_data, 1, min)
pdt_data = pdt_data / apply(pdt_data, 1, max)

beautiful_result_norm = reshape2::melt(data=pdt_data)

colnames(beautiful_result_norm) = c("GeneNames", "Pseudotime", "ExpressionValue")

print("preparing to plot genes by expression level")

plot.genes = ggplot(data = beautiful_result_norm, aes(x = Pseudotime, y = GeneNames))
plot.genes = plot.genes + geom_tile(aes(fill = ExpressionValue),  width=1.001, height=1.001)
plot.genes = plot.genes + scale_fill_gradient2(low = "deepskyblue", high = "firebrick3", mid = "darkolivegreen3", midpoint = 0.5, name = "Minmax normalized gene expression")
plot.genes = plot.genes + theme(legend.position = "bottom", legend.text = element_text(size = 25, angle = 90),
                                 legend.title = element_text(size = 25),
                                 legend.key.width = unit(2, "cm"),
                                 axis.text.x = element_blank(), axis.title.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 8))

pdf(file.path(output_folder, "expression_vs_norm_expression.pdf"), width = 13, height = 35)
plot_grid(plot.density, plot.genes, nrow = 2, align = "v", rel_heights = c(1/9, 8/9))
dev.off()

print("plotted expression_vs_norm_expression.pdf in output folder")

# plot non-normalized expression
###################################################################################################
raw_data_genes = as.matrix(seurat.obj@data[rev(genes.in.order), order(pseudotime$Pseudotime)])
print("made raw_data_genes matrix")
raw_data_genes = t(apply(raw_data_genes, 1, adaptive.moving_average, kernel = 15, minim_kernel = 1, range.factor=15))
print("raw_data_genes matrix has applied apply function and t")

# group by pdt
pdt = range(pseudotime$Pseudotime)
pdt = seq(pdt[1], pdt[2], length.out=100)
pdt_data = c()
for (k in 1:nrow(raw_data_genes)){
  for(j in 1:length(pdt)){
    local_pdt = pdt[j]
    pdt_index = abs(pseudotime$Pseudotime[order(pseudotime$Pseudotime)] - local_pdt)
    pdt_index = which(pdt_index == min(pdt_index))
    pdt_data = c(pdt_data, raw_data_genes[k, pdt_index])
  }
}
pdt_data = matrix(data=pdt_data, nrow=nrow(raw_data_genes), byrow=T)
rownames(pdt_data) = rownames(raw_data_genes)
#colnames(pdt_data) = paste("PDT", 1:ncol(pdt_data), sep = "")
colnames(pdt_data) = paste("PDT", 1:100, sep = "")

# smooth a bit the pdt_data matrx
pdt_data = t(apply(pdt_data, 1, ma, kernel = 7))

beautiful_result_nonnorm = reshape2::melt(data=pdt_data)
colnames(beautiful_result_nonnorm) = c("GeneNames", "Pseudotime", "ExpressionValue")

plot.genes = ggplot(data = beautiful_result_nonnorm, aes(x = Pseudotime, y = GeneNames))
plot.genes = plot.genes + geom_tile(aes(fill = ExpressionValue),  width=1.001, height=1.001)
plot.genes = plot.genes + scale_fill_gradient2(low = "deepskyblue", high = "firebrick3", mid = "darkolivegreen3", midpoint = mean(range(pdt_data)), name = "Gene expression")
plot.genes = plot.genes + theme(legend.position = "bottom", legend.text = element_text(size = 25, angle = 90),
                                 legend.title = element_text(size = 25),
                                 legend.key.width = unit(2, "cm"),
                                 axis.text.x = element_blank(), axis.title.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 8))

pdf(file.path(output_folder, "expression_vs_nonnorm_expression.pdf"), width = 13, height = 35)
plot_grid(plot.density, plot.genes, nrow = 2, align = "v", rel_heights = c(1/9, 8/9))
dev.off()

print("plotted genes by expression_vs_nonnorm_expression.pdf")

# save diffusion map coordinates and expression data for found genes
by.pdt.order = order(pseudotime$Pseudotime)

dm.df = read.csv(file.path(output_folder_material, "dm.csv"), row.names = 1, header = F)
dm.df = as.data.frame(dm.df[, 1:3])
dm.df$Labels = factor(seurat.obj@ident, levels = cell.labels)
dm.df$Colours = mapvalues(x = dm.df$Labels, from = cell.labels, to = cell.colours)
dm.df = dm.df[by.pdt.order, ]
colnames(dm.df) = c("DM1", "DM2", "DM3", "Labels", "Colours")

print("writing pdt_and_expression.csv")

expression_data_and_pdt = as.data.frame(t(as.matrix(seurat.obj@data[rev(genes.in.order), by.pdt.order])))
pdt.data = data.frame(Pseudotime = pseudotime[by.pdt.order, c(1)])
pdt.data = cbind(dm.df, pdt.data, expression_data_and_pdt)
pdt.data.fp = file.path(output_folder, "pdt_and_expression.csv")
write.csv(pdt.data, pdt.data.fp, row.names = F)

# make interactive diffusion map
command = sprintf("%s html_3D_viewer_and_plotter.py %s %s", python.addr, file.path(output_folder, "Interactive_Pseudotime.html"), pdt.data.fp)
system(command, wait = T)

# save the plotting material, just in case
plot.data.objects = list(pseudotime = pseudotime, beautiful_result_norm = beautiful_result_norm, beautiful_result_nonnorm = beautiful_result_nonnorm)
saveRDS(plot.data.objects, file.path(output_folder, "ploting_material.RDS"))

unlink(output_folder_material, recursive=T, force=T)

print("Ended beautifully ... ")
