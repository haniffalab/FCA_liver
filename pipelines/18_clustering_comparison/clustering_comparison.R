args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg    = args[1]
set.ident.arg      = args[2]
n_clusters.arg     = args[3]
type.to.colour.arg = args[4]
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

output_folder_material = file.path(output_folder, "material")
dir.create(output_folder_material)

seurat.addr = file.path("../../data", seurat.addr)

n_clusters = as.numeric(n_clusters)

source("../../tools/bunddle_utils.R")

library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

seurat.obj = SetAllIdent(seurat.obj, id=set.ident)

# saving pca data to disk
pca.df = seurat.obj@dr$pca@cell.embeddings
pca.df = cbind(pca.df, data.frame(Ident = as.vector(seurat.obj@ident)))
pca.fname = file.path(output_folder_material, 'pca.csv')
write.csv(pca.df, pca.fname)

print("Performing clustering ... ")
command = sprintf("%s ./clustering.py %s %s", python.addr, output_folder, n_clusters)
system(command, wait = T)
agg.clusters = read.csv(file.path(output_folder_material, 'agglomerative_clustering.csv'))
colnames(agg.clusters) = c('CellNames', 'ClusterIdx')
gau.clusters = read.csv(file.path(output_folder_material, 'gaussian_mixture.csv'))
colnames(gau.clusters) = c('CellNames', 'ClusterIdx')
metrics_file = file.path(output_folder_material, 'agreement_measures.txt')
metrics_file = file(metrics_file)
metrics = readLines(metrics_file)
close(metrics_file)
metrics = round(as.numeric(metrics), digits=2)
metrics = matrix(metrics, nrow=2)
colnames(metrics) = c('Rand Index', 'Adj. Mutual Information')
rownames(metrics) = c('Louvain vs Agglomerative', 'Louvain vs Gaussian Mixture')
metrics = as.table(metrics)
library(gridExtra)
library(grid)
grid.fname = file.path(output_folder, 'metrics.pdf')
pdf(grid.fname, width = 6, height = 1.5)
grid.table(metrics)
dev.off()

seurat.obj@meta.data$gmm = gau.clusters$ClusterIdx
seurat.obj@meta.data$agg = agg.clusters$ClusterIdx


if (!is.na(type.to.colour)){
  type.to.colour = read.csv(file.path("../../resources", type.to.colour))
  filter.key   =  type.to.colour$CellTypes %in% as.vector(unique(seurat.obj@ident))
  cell.labels  = as.vector(type.to.colour$CellTypes[filter.key])
  cell.colours = as.vector(type.to.colour$Colours[filter.key])
}else{
  cell.labels  = sort(as.vector(unique(seurat.obj@ident)))
  cell.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
}

louvain_idx = as.vector(seurat.obj@ident)
gmm_idx = as.vector(seurat.obj@meta.data$gmm)
agg_idx = as.vector(seurat.obj@meta.data$agg)
umap_data = as.data.frame(seurat.obj@dr$umap@cell.embeddings)

plot.umap.louvain = dr.plot(point.labels=louvain_idx, dr1=umap_data$umap1, dr2=umap_data$umap2, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, use.labels=cell.labels)
plot.umap.louvain = plot.umap.louvain + ggtitle('Cell types from\nLouvain clustering')
gmm_types = as.vector(unique(gmm_idx))
gmm_types = cell.labels[cell.labels %in% gmm_types]
gmm_colours = mapvalues(x=gmm_types, from=cell.labels, to=cell.colours)
gmm_index   = seq(1, length(cell.labels))[cell.labels %in% gmm_types]
plot.umap.gm = dr.plot(point.labels=gmm_idx, dr1=umap_data$umap1, dr2=umap_data$umap2, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=gmm_colours, use.labels=gmm_types, index.map=gmm_index)
plot.umap.gm = plot.umap.gm + ggtitle('Cell types from\nGaussian Mixture clustering')
agg_types = as.vector(unique(agg_idx))
agg_types = cell.labels[cell.labels %in% agg_types]
agg_colours = mapvalues(x=agg_types, from=cell.labels, to=cell.colours)
agg_index   = seq(1, length(cell.labels))[cell.labels %in% agg_types]
plot.umap.agg = dr.plot(point.labels=agg_idx, dr1=umap_data$umap1, dr2=umap_data$umap2, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=agg_colours, use.labels=agg_types, index.map=agg_index)
plot.umap.agg = plot.umap.agg + ggtitle('Cell types from\nAgglomerative clustering')

dm.fname = file.path(output_folder, 'clustering_comparisons.pdf')
pdf(dm.fname, width = 18, height = 5)
plot_grid(plot.umap.louvain, plot.umap.gm, plot.umap.agg, nrow=1)
dev.off()

n.cols = min(2, length(cell.labels))
n.rows = ceiling(length(cell.labels) / n.cols)
plot.legend = plot.indexed.legend(label.vector=cell.labels, color.vector=cell.colours, ncols=n.cols, left.limit=.2, symbol.size=10, text.size=6, padH=.6, padV=.6)
legend.file.name = file.path(output_folder, 'legend.pdf')
pdf(legend.file.name, width = 1.5 + .15 * n.cols * max(unlist(lapply(cell.labels, nchar))), height = .5 + n.rows * .35)
print(plot.legend)
dev.off()

plot.umap.louvain = dr.plot(point.labels=louvain_idx, dr1=umap_data$umap1, dr2=umap_data$umap2, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, use.labels=cell.labels, annotate.plot=F)
plot.umap.louvain = plot.umap.louvain + ggtitle('Cell types from\nLouvain clustering')
gmm_types = as.vector(unique(gmm_idx))
gmm_types = cell.labels[cell.labels %in% gmm_types]
gmm_colours = mapvalues(x=gmm_types, from=cell.labels, to=cell.colours)
gmm_index   = seq(1, length(cell.labels))[cell.labels %in% gmm_types]
plot.umap.gm = dr.plot(point.labels=gmm_idx, dr1=umap_data$umap1, dr2=umap_data$umap2, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=gmm_colours, use.labels=gmm_types, index.map=gmm_index, annotate.plot=F)
plot.umap.gm = plot.umap.gm + ggtitle('Cell types from\nGaussian Mixture clustering')
agg_types = as.vector(unique(agg_idx))
agg_types = cell.labels[cell.labels %in% agg_types]
agg_colours = mapvalues(x=agg_types, from=cell.labels, to=cell.colours)
agg_index   = seq(1, length(cell.labels))[cell.labels %in% agg_types]
plot.umap.agg = dr.plot(point.labels=agg_idx, dr1=umap_data$umap1, dr2=umap_data$umap2, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=agg_colours, use.labels=agg_types, index.map=agg_index, annotate.plot=F)
plot.umap.agg = plot.umap.agg + ggtitle('Cell types from\nAgglomerative clustering')

dm.fname = file.path(output_folder, 'clustering_comparisons_no_dots.pdf')
pdf(dm.fname, width = 18, height = 5)
plot_grid(plot.umap.louvain, plot.umap.gm, plot.umap.agg, nrow=1)
dev.off()

unlink(output_folder_material, recursive=T, force=T)

print("Ended beautifully ... ")
