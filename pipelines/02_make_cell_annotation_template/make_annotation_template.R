args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg    = args[1]
clustering.res.arg = args[2]
DE.downsample.arg  = args[3]
sample.arg         = args[4]
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

source("../../tools/bunddle_utils.R")

library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(BiocParallel)

dr.plot.indexed.clusters <- function(point.labels, dr1, dr2, dr1.name, dr2.name, no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2){
  df.dr <- data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2)
  p.labels <- unique(as.vector(point.labels))
  p.labels <- as.character(sort(as.numeric(p.labels)))
  p.labels.medians <- aggregate(df.dr[, 2:3], list(df.dr$Cell.Labels), median)
  set.seed(random_state)
  plt.colours <- sample(colorRampPalette(brewer.pal(12, "Paired"))(length(p.labels)))
  index.map <- p.labels
  plot.obj <- ggplot(data = df.dr, aes(x = DR1, y = DR2, color = Cell.Labels))
  plot.obj <- plot.obj + geom_point(size = pt.size)
  plot.obj <- plot.obj + scale_color_manual(values=plt.colours)
  plot.obj <- plot.obj + stat_density2d(geom="density2d", aes(x=DR1, y=DR2,alpha=5/10), size=.2, contour=TRUE,bins=7,h=1.5)
  plot.obj <- plot.obj + geom_point(data=p.labels.medians,aes(x = DR1, y = DR2), colour = "gray", size = plt.lb.sz, fill = plt.colours, alpha = .5, pch = 21)
  plot.obj <- plot.obj + annotate("text", x=p.labels.medians$DR1, y = p.labels.medians$DR2, label = as.vector(p.labels.medians$Group.1), size = txt.lb.size)
  if (no.legend){
    plot.obj <- plot.obj + theme(legend.position="none")
  }else{
    plot.obj <- plot.obj + guides(color = guide_legend(override.aes = list(size=5)))
  }
  plot.obj <- plot.obj + xlab(dr1.name) + ylab(dr2.name)
  return(plot.obj)
}

dr.plot <- function(point.labels, dr1, dr2, dr1.name, dr2.name, no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2, use.cols = FALSE, index.map = c()){
  df.dr <- data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2)
  p.labels <- sort(unique(as.vector(point.labels)))
  df.dr$Cell.Labels <- factor(df.dr$Cell.Labels, levels=p.labels)
  p.labels.medians <- aggregate(df.dr[, 2:3], list(df.dr$Cell.Labels), median)
  df.dr$Cell.Labels <- mapvalues(x = df.dr$Cell.Labels, from = p.labels, to = paste(1:length(p.labels), p.labels, sep = " "))
  if(is.logical(use.cols)){
    set.seed(random_state)
    plt.colours <- sample(colorRampPalette(brewer.pal(12, "Paired"))(length(p.labels)))
    index.map <- 1:length(p.labels)
  }else{
    plt.colours <- use.cols
  }
  plot.obj <- ggplot(data = df.dr, aes(x = DR1, y = DR2, color = Cell.Labels))
  plot.obj <- plot.obj + geom_point(size = pt.size)
  plot.obj <- plot.obj + scale_color_manual(values=plt.colours)
  #plot.obj <- plot.obj + stat_density2d(geom="density2d", aes(x=DR1, y=DR2,alpha=5/10), size=.2, contour=TRUE,bins=7,h=1.5)
  plot.obj <- plot.obj + geom_point(data=p.labels.medians,aes(x = DR1, y = DR2), colour = "gray", size = plt.lb.sz, fill = "gray", alpha = .5, pch = 21)
  plot.obj <- plot.obj + annotate("text", x=p.labels.medians$DR1, y = p.labels.medians$DR2, label = index.map, size = txt.lb.size)
  if (no.legend){
    plot.obj <- plot.obj + theme(legend.position="none")
  }else{
    plot.obj <- plot.obj + guides(color = guide_legend(override.aes = list(size=5)))
  }
  plot.obj <- plot.obj + xlab(dr1.name) + ylab(dr2.name)
  return(plot.obj)
}

dr.plot.group <- function(point.labels, dr1, dr2, dr1.name, dr2.name, group.name, pt.size = .4){
  df.dr <- data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2)
  p.labels <- sort(unique(as.vector(point.labels)))
  df.dr$Cell.Labels <- factor(df.dr$Cell.Labels, levels=p.labels)
  group.index <- which(p.labels == group.name)
  plt.colours <- rep("#bae1ff", length(p.labels))
  plt.colours[group.index] <- "#0D7D75"
  plot.obj <- ggplot(data = df.dr, aes(x = DR1, y = DR2, color = Cell.Labels))
  plot.obj <- plot.obj + geom_point(size = pt.size)
  plot.obj <- plot.obj + scale_color_manual(values=plt.colours)
  plot.obj <- plot.obj + theme(legend.position="none")
  plot.obj <- plot.obj + xlab(dr1.name) + ylab(dr2.name) + ggtitle(group.name)
  return(plot.obj)
}

tabulate.seurat.by.cluster <- function(seurat.obj, slot1, slot2, save.at, width, height, saveas.pdf = F){
  "used to build tables that show contingency distribution of cells by 2 different labeling criteria"
  "these are slot1 and slot2 which should be in the meta.data slot of the seurat object"
  for (i in 1:length(levels(seurat.obj@ident))){
    cluster = levels(seurat.obj@ident)[i]
    cells.cluster <- colnames(seurat.obj@data)[seurat.obj@ident == cluster]
    cells.indices <- match(cells.cluster, colnames(seurat.obj@data))
    base.com <- paste(substitute(seurat.obj), "meta.data", sep = "@")
    com1 <- paste(base.com, slot1, sep = "$")
    com1 <- sprintf("%s[cells.indices]", com1)
    com2 <- paste(base.com, slot2, sep = "$")
    com2 <- sprintf("%s[cells.indices]", com2)
    command <- sprintf("tally <- table(%s, %s)", com1, com2)
    eval(parse(text = command))
    command <- sprintf("tally.rez <- cbind(tally, `Total by %s` = rowSums(tally))", slot1)
    eval(parse(text = command))
    command <- sprintf("tally.rez <- rbind(tally.rez, `Total by %s` = c(colSums(tally), length(cells.cluster)))", slot2)
    eval(parse(text = command))
    print(tally.rez)
    if (saveas.pdf){
      filename <- paste(paste("tally_", cluster, sep = ""), ".pdf", sep = "")
      filename <- file.path(save.at, filename)
      pdf(filename, width = width, height = height)
      grid.table(tally.rez)
      dev.off()
    }else{
      filename <- paste(paste("tally_", cluster, sep = ""), ".png", sep = "")
      filename <- file.path(save.at, filename)
      png(filename, width = width, height = height)
      grid.table(tally.rez)
      dev.off()
    }
  }
}

FindMarker.wrapper <- function(markers.for){
  if (DE.downsample){
    markers = FindMarkers(seurat.obj_d, ident.1=markers.for, only.pos = F, min.pct=0.25, genes.use=rownames(seurat.obj_d@data),
                          thresh.use = 0.25, test.use = "wilcox", random.seed = 42, print.bar=T, do.print=T)
  }else{
    markers = FindMarkers(seurat.obj, ident.1=markers.for, only.pos = F, min.pct=0.25, genes.use=rownames(seurat.obj@data),
                          thresh.use = 0.25, test.use = "wilcox", random.seed = 42, print.bar=T, do.print=T)
  }
  markers$cluster = markers.for
  markers$gene = rownames(markers)
  markers
}

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

# process the data (normalize - scale - variable genes - pca - tsne)
print("Normalizing data ... ")
seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
print("Computing variable genes ... ")

# find all clusters
print("Clustering data ... ")
seurat.obj <- FindClusters(object = seurat.obj, reduction.type = "pca",
                           dims.use = 1:20, resolution = clustering.res, save.SNN = T, algorithm=1)
print(paste("Number of clusters: ", print(length(levels(seurat.obj@ident))), sep = ""))
seurat.obj@meta.data$LouvainClustering = as.vector(seurat.obj@ident)
print(table(seurat.obj@meta.data$LouvainClustering))

print("Saving Seurat object")
saveRDS(seurat.obj, seurat.addr)
print('Seurat object saved')

# writing marker genes to disk
if (DE.downsample){
  cluster.ids <-unique(as.vector(seurat.obj@ident))
  cells.to.keep <- c()
  for (k in 1:length(cluster.ids)){
    cluster.id <- cluster.ids[k]
    cell.ids <- names(seurat.obj@ident)[seurat.obj@ident == cluster.id]
    cell.ids <- which(names(seurat.obj@ident) %in% cell.ids )
    cells.to.keep <- c(sample(x=cell.ids, size=min(300, length(cell.ids)), replace=F), cells.to.keep)
  }
  seurat.obj_d <- SubsetData(object=seurat.obj, cells.use=names(seurat.obj@ident)[cells.to.keep])
  seurat.obj_d <- NormalizeData(object = seurat.obj_d, normalization.method = "LogNormalize", scale.factor = 10000)
  print("Calculating marker genes: finished subseting, currently actually calculating the markers ... ")
  
  Markers <- bplapply(sort(as.vector(unique(seurat.obj_d@meta.data$LouvainClustering))), FindMarker.wrapper,BPPARAM=MulticoreParam(5))
}else{
  print("Calculating marker genes: finished subseting, currently actually calculating the markers ... ")
  Markers <- bplapply(sort(as.vector(unique(seurat.obj@meta.data$LouvainClustering))), FindMarker.wrapper,BPPARAM=MulticoreParam(5))
}
marker.genes = Reduce(f=rbind, x=Markers)

print("Saving marker genes ... ")
write.csv(marker.genes, file.path(output_folder, "all_markers.csv"))

print('Creating and saving to disk annotation marker genes')

gene_db = read.csv('./gene_info.csv')
rownames(gene_db) = as.vector(gene_db$gene.symbol)
marker.genes.top = marker.genes %>% group_by(cluster) %>% top_n(50, avg_logFC)
gene_to_pop = read.csv("./gene_to_pop.tsv", sep = '\t', header = F)
colnames(gene_to_pop) = c('Gene', 'Population')

marker.genes.unique = unique(as.vector(marker.genes.top$gene))
gene_info = gene_db[marker.genes.unique, ]

# get gene name
gene.name = mapvalues(x=as.vector(marker.genes.top$gene), from=as.vector(gene_info$gene.symbol), 
                      to=as.vector(gene_info$gene.name))
marker.genes.top = cbind(as.data.frame(marker.genes.top), data.frame(GeneName = gene.name))

# get also present in 
also_present_in = function(gene.sym){
  part = as.vector(marker.genes.top[as.vector(marker.genes.top$gene) == gene.sym, ]$cluster)
  if (length(part) > 1){
    return(paste(part, collapse = ', '))
  }else{
    return('')
  }
}
present_in = unlist(lapply(as.list(marker.genes.unique), also_present_in))
present_in = mapvalues(x=as.vector(marker.genes.top$gene), from=as.vector(marker.genes.unique), to=as.vector(present_in))
marker.genes.top = cbind(marker.genes.top, data.frame(AlsoPresentInClusters = present_in))

# get cell type flag
cell_type_flag = c()
for(i in 1:dim(marker.genes.top)[1]){
  gene.sym = as.vector(marker.genes.top$gene)[i]
  pops = as.vector(gene_to_pop$Population)[as.vector(gene_to_pop$Gene) == gene.sym]
  if(length(pops) == 0){
    pops = ''
  }
  cell_type_flag = c(cell_type_flag, pops)
}
marker.genes.top = cbind(marker.genes.top, data.frame(CellTypeFlag = cell_type_flag))

# get gene summary
gene.summary = mapvalues(x=as.vector(marker.genes.top$gene), from=as.vector(gene_info$gene.symbol), 
                         to=as.vector(gene_info$gene.summary))
marker.genes.top = cbind(as.data.frame(marker.genes.top), data.frame(Summary = gene.summary))

# get reactom
reactome.pathway = mapvalues(x=as.vector(marker.genes.top$gene), from=as.vector(gene_info$gene.symbol), 
                             to=as.vector(gene_info$reactome.pathway))
marker.genes.top = cbind(as.data.frame(marker.genes.top), data.frame(Reactom = reactome.pathway))

# get gene family
gene.family = mapvalues(x=as.vector(marker.genes.top$gene), from=as.vector(gene_info$gene.symbol), 
                        to=as.vector(gene_info$gene.family))
marker.genes.top = cbind(as.data.frame(marker.genes.top), data.frame(GeneFamily = gene.family))

write.csv(marker.genes.top, file.path(output_folder, "annotation_markers.csv"))

update.template <- data.frame(Cluster = sort(as.vector(unique(seurat.obj@ident))), Identity = rep("None", length(unique(seurat.obj@ident))))
write.csv(update.template, file.path(output_folder, "update_template.csv"), row.names = F)

print("compiling the template")

df <- data.frame(CellNames = names(seurat.obj@ident),
                 ClusterIndex = as.vector(seurat.obj@ident),
                 tSNEx = seurat.obj@dr$tsne@cell.embeddings[, 1],
                 tSNEy = seurat.obj@dr$tsne@cell.embeddings[, 2],
                 UMAPx = seurat.obj@dr$umap@cell.embeddings[, 1],
                 UMAPy = seurat.obj@dr$umap@cell.embeddings[, 2],
                 Sample = seurat.obj@meta.data[, sample])

# transfer the compile.py file to the output
file.copy(from='compile_template.py', to=file.path(output_folder, 'compile_template.py'))

CURDIR = getwd()
setwd(output_folder)

dir.create("graphs")

# make tsne and umap plots by clusters
plot.tsne <- dr.plot.indexed.clusters(point.labels=df$ClusterIndex, dr1=df$tSNEx, dr2=df$tSNEy, dr1.name="tSNE-x", dr2.name="tSNE-y", no.legend = T, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2)
plot.umap <- dr.plot.indexed.clusters(point.labels=df$ClusterIndex, dr1=df$UMAPx, dr2=df$UMAPy, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend = T, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2)
png("./graphs/dr.png", width = 1200, height = 700)
plot_grid(plot.tsne, plot.umap)
dev.off()

# make tsne and umap plots by sample
plot.tsne <- dr.plot(point.labels=df$Sample, dr1=df$tSNEx, dr2=df$tSNEy, dr1.name="tSNE-x", dr2.name="tSNE-y", no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = 1)
plot.umap <- dr.plot(point.labels=df$Sample, dr1=df$UMAPx, dr2=df$UMAPy, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = 1)
png("./graphs/dr_sample.png", width = 1200, height = 700)
plot_grid(plot.tsne, plot.umap)
dev.off()

# create cell tally plots, tsne plots and umap plots
no.clusters <- length(levels(seurat.obj@ident))
for (i in 1:no.clusters){
  cluster.name <- levels(seurat.obj@ident)[i]
  plot.tsne <- dr.plot.group(point.labels=df$ClusterIndex, dr1=df$tSNEx, dr2=df$tSNEy, dr1.name="tSNE1", dr2.name="tSNE2", group.name=cluster.name, pt.size = .4)
  plot.umap <- dr.plot.group(point.labels=df$ClusterIndex, dr1=df$UMAPx, dr2=df$UMAPy, dr1.name="UMAP1", dr2.name="UMAP2", group.name=cluster.name, pt.size = .4)
  graph.addr <- paste(paste("cluster_dr_", cluster.name, sep = ""), ".png", sep = "")
  graph.addr <- file.path("graphs", graph.addr)
  png(graph.addr, width = 400, height = 900)
  print(plot_grid(plot.tsne, plot.umap, nrow = 2))
  dev.off()
  print(cluster.name)
}

# plot cell numbers by sample and gate for all the clusters ("sort.ids", "fetal.ids")
tabulate.seurat.by.cluster(seurat.obj, "tissue", "tissue", save.at="./graphs", width=1110, height=110, saveas.pdf = F)

# compile template annotation powerpoint
system(paste(python.addr, "compile_template.py", sep = " "), wait = T)

setwd(CURDIR)

print('Finished.')
