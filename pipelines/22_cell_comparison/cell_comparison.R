args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.query.addr.arg = args[1]
seurat.ref.addr.arg   = args[2]
set.ident.query.arg   = args[3]
set.ident.ref.arg     = args[4]
cell.types.query.arg  = args[5]
cell.types.ref.arg    = args[6]
dims.plot.arg         = args[7]
compute.DEGs.arg      = args[8]
"

DEG.number.limit = 20

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
label1 = gsub(pattern="\\.RDS", replacement="", x=seurat.query.addr)
label2 = gsub(pattern="\\.RDS", replacement="", x=seurat.ref.addr)
output_folder = paste(output_folder, label1, label2, sep = "_")
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

seurat.query.addr = file.path("../../data", seurat.query.addr)
seurat.ref.addr   = file.path("../../data", seurat.ref.addr)

source("../../tools/bunddle_utils.R")

gene_mean_expression = function(seurat.obj){
  no.genes = nrow(seurat.obj@data)
  start_index = 1
  while (start_index < no.genes){
    end_index = start_index + 999
    end_index = min(end_index, no.genes)
    expression.data_ = data.matrix(seurat.obj@data[start_index:end_index, ])
    expression.data_ = t(expression.data_)
    expression.data_ = as.data.frame(expression.data_)
    expression.data_ = cbind(data.frame(CellLabels = as.vector(seurat.obj@ident)), expression.data_)
    expression.data_ = aggregate(expression.data_[2:dim(expression.data_)[2]], list(expression.data_$CellLabels), mean)
    expression.data_ = cbind(data.frame(CellType = expression.data_$Group.1), expression.data_[, 2:dim(expression.data_)[2]])
    rownames(expression.data_) = expression.data_$CellType
    expression.data_ = expression.data_[, 2:ncol(expression.data_)]
    print(start_index)
    if (start_index == 1){
      expression.data = expression.data_
    }else{
      expression.data = cbind(expression.data, expression.data_)
    }
    start_index = start_index + 1000
  }
  expression.data
}

make_heatmap = function(table.data, measure.name, low.col = "darkslateblue", mid.col = "floralwhite", high.col = "goldenrod1"){
  expression.melt = reshape::melt(data=as.matrix(table.data))
  colnames(expression.melt) = c("CellTypesQuery", "CellTypesReference", measure.name)
  expression.melt$CellTypesQuery     = factor(as.vector(expression.melt$CellTypesQuery),     levels=rev(query.labels))
  expression.melt$CellTypesReference = factor(as.vector(expression.melt$CellTypesReference), levels = ref.labels)
  heatmap.plot = ggplot(expression.melt, aes(CellTypesReference, CellTypesQuery, fill = expression.melt[, 3])) +
    geom_tile(color = "black") +
    scale_fill_gradientn(colours = c(low.col, mid.col, high.col), values = c(0, .5, 1.0)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                                      axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill=measure.name)
  heatmap.plot
}

plot_DEGs = function(query.pop){
  print("running plot_degs function")
  vln.grobs = c()
  first_time = T
  for (i in seq_along(ref.labels)){
    ref.pop = ref.labels[i]
    # compute DEGs for query vs refs
    if (query.pop != ref.pop){
      print("printing query pop")
      print(query.pop)
      print("printing ref.pop")
      print(ref.pop)
      #DEGs = FindMarkers(object=seurat.obj, ident.1=query.pop, ident.2=ref.pop, genes.use=rownames(seurat.obj@raw.data), only.pos=F, max.cells.per.ident = 500)
      DEGs = FindMarkers(object=seurat.obj, ident.1=query.pop, ident.2=ref.pop, features=rownames(seurat.obj@raw.data), only.pos=F, max.cells.per.ident = 500)
      print("printing DEGs that are found straight after findmarkers is run")
      print(DEGs)
      DEGs$Gene = rownames(DEGs)
      DEGs$Reference_Population = gsub(pattern="R::", replacement="", x=ref.pop)
      DEGs$Query_Population     = gsub(pattern="Q::", replacement="", x=query.pop)
      if (first_time){
        DEG.matrix = DEGs
        first_time = F
      }else{
        DEG.matrix = rbind(DEG.matrix, DEGs)
      }
      
      print("printing DEGs after if loop to add cols")
      print(DEGs)
      pos.markers = DEGs %>% top_n(deg.span, avg_logFC)
      pos.markers = as.vector(pos.markers$Gene)
      neg.markers = DEGs %>% top_n(-deg.span, avg_logFC)
      neg.markers = rev(as.vector(neg.markers$Gene))
      markers     = c(pos.markers, neg.markers)
      
      # make marker data - used downstream for plotting
      markers.data.query = t(as.matrix(seurat.obj@data[markers, as.vector(seurat.obj@ident) == query.pop]))
      markers.data.query = reshape::melt(markers.data.query)
      colnames(markers.data.query) = c("Population", "Gene", "Expression")
      markers.data.query$Population = "Query_Population"
      
      markers.data.ref   = t(as.matrix(seurat.obj@data[markers, as.vector(seurat.obj@ident) == ref.pop]))
      markers.data.ref = reshape::melt(markers.data.ref)
      colnames(markers.data.ref) = c("Population", "Gene", "Expression")
      markers.data.ref$Population = "Reference_Population"
      
      markers.data = as.data.frame(rbind(markers.data.query, markers.data.ref))
      markers.data$Gene = factor(as.vector(markers.data$Gene), levels = markers)
            
      # make violin plot object
      vln.plot = ggplot(markers.data, aes(x = Population, y = Expression, col = Population, fill = Population))
      vln.plot = vln.plot + geom_jitter(size = .3)
      vln.plot = vln.plot + facet_wrap(~Gene, nrow = 1)
      vln.plot = vln.plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position="none")
      vln.plot = vln.plot + theme(axis.title.y = element_text(angle = 90, hjust = 1)) + ylab(ref.pop)
      
      unique_vln_plot = paste("VlnPlot", i, sep = "_")
      eval(parse(text = sprintf('%s = vln.plot', unique_vln_plot)))
      vln.grobs = c(vln.grobs, unique_vln_plot)
    }
  }
  vln.grobs = paste(vln.grobs, collapse = ", ")
  eval(parse(text = sprintf('vln.grobs = list(%s)', vln.grobs)))
  
  DEG.fname = gsub(pattern="/", replacement='_', x=query.pop, fixed = T)
  DEG.fname = gsub(pattern="Q::", replacement="", x=DEG.fname)
  DEG.fname = paste(DEG.fname, ".png", sep = "")
  DEG.fname = file.path(deg.folder, DEG.fname)
  
  print(DEG.fname)
  
  png(DEG.fname, width = plot.W, height = plot.H)
  grid.arrange(grobs = vln.grobs, ncol = 1)
  dev.off()
  
  DEG.fname = gsub(pattern="/", replacement='_', x=query.pop, fixed = T)
  DEG.fname = gsub(pattern="Q::", replacement="", x=DEG.fname)
  DEG.fname = paste(DEG.fname, ".csv", sep = "")
  DEG.fname = file.path(deg.folder, DEG.fname)
  
  print(DEG.fname)
  
  write.csv(DEG.matrix, DEG.fname)
}

DownSample = function(seurat.obj, limit){
  pop.types = as.vector(unique(seurat.obj@ident))
  to.keep = c()
  for(k in seq_along(pop.types)){
    pop.type = pop.types[k]
    pop.names = names(seurat.obj@ident)[as.vector(seurat.obj@ident) == pop.type]
    if (length(pop.names) > limit){
      pop.names = sample(x=pop.names, size=limit, replace=FALSE)
    }
    to.keep = c(to.keep, pop.names)
  }
  SubsetData(object=seurat.obj, cells.use=to.keep)
}

library(Seurat)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(gridExtra)
library(BiocParallel)

#######################################################################################################

# load data
print("loading query data set ... ")
seurat.query.obj = readRDS(seurat.query.addr)
seurat.query.obj = SetAllIdent(object=seurat.query.obj, id=set.ident.query)
if (ncol(seurat.query.obj@data) > 100000){
  print("Query data too big. Downsampling query data ...")
  seurat.query.obj = DownSample(seurat.query.obj, limit = 2000)
}
print("Query data loaded.")

if(seurat.ref.addr != seurat.query.addr){
  print("loading reference data set ... ")
  seurat.ref.obj = readRDS(seurat.ref.addr)
  seurat.ref.obj   = SetAllIdent(object=seurat.ref.obj,   id=set.ident.ref)
  if(ncol(seurat.ref.obj@data) > 100000){
    print("Reference data too big. Downsampling reference data ...")
    seurat.ref.obj = DownSample(seurat.ref.obj, limit = 2000)
  }
  print("Reference data loaded.")
}else{
  print("Query and reference data sets are the same. Copying the query to make the reference data sets ...")
  if (set.ident.query != set.ident.ref){
    stop("The reference and query data sets are the same but the set identities are different. Terminating script ...")
  }
  seurat.ref.obj = seurat.query.obj
  seurat.obj     = seurat.query.obj
}

# load lists of cell types
if (is.na(cell.types.query)){
  cell.types.query = sort(as.vector(unique(seurat.query.obj@ident)))
}else{
  # open file connection
  cell.types.query.file = file.path("../../resources", cell.types.query)
  cell.types.query.file = file(cell.types.query.file, "r")
  cell.types.query = readLines(cell.types.query.file)
  close(cell.types.query.file)
}

if(is.na(cell.types.ref)){
  cell.types.ref = sort(as.vector(unique(seurat.ref.obj@ident)))
}else{
  # open file connection
  cell.types.ref.file = file.path("../../resources/", cell.types.ref)
  cell.types.ref.file = file(cell.types.ref.file, "r")
  cell.types.ref = readLines(cell.types.ref.file)
  close(cell.types.ref.file)
}

# subset seurat objects
if (!all(as.vector(unique(seurat.query.obj@ident)) %in% cell.types.query)){
  print("Subsetting the query data ...")
  print("as.vector(unique(seurat.query.obj@ident)")
  print(as.vector(unique(seurat.query.obj@ident)))
  print("cell types query")
  print(cell.types.query)
  seurat.query.obj = SubsetData(object=seurat.query.obj, ident.use=cell.types.query, do.clean=T, subset.raw=T)
}
if(!all(as.vector(unique(seurat.ref.obj@ident)) %in% cell.types.ref)){
  print("Subsetting the reference data ...")
  seurat.ref.obj   = SubsetData(object=seurat.ref.obj,   ident.use=cell.types.ref,   do.clean=T, subset.raw=T)
}

#remove any cell type represented by to few cells to allow a statistical test of DEGs
if (compute.DEGs){
  
cell.types.query.fail = names(table(seurat.query.obj@ident)[table(seurat.query.obj@ident) < DEG.number.limit])
cell.types.ref.fail   = names(table(seurat.ref.obj@ident)[table(seurat.ref.obj@ident) < DEG.number.limit])

if (length(cell.types.query.fail) > 0){
  print(sprintf("The query data has cell population represented by too few cells to allow DEG statistical test. These are: '%s';", paste(cell.types.query.fail, collapse = ",")))
  print("Removing these cell types from query data ...")
  cell.types.query = cell.types.query[!(cell.types.query %in% cell.types.query.fail)]
  seurat.query.obj = SubsetData(object=seurat.query.obj, ident.use=cell.types.query, do.clean=T, subset.raw=T)
}

if (length(cell.types.ref.fail) > 0){
  print(sprintf("The reference data has cell population represented by too few cells to allow DEG statistical test. These are:' %s'; And will be removed.", paste(cell.types.ref.fail, collapse = ",")))
  print("Removing these cell types from reference data ...")
  cell.types.ref = cell.types.ref[!(cell.types.ref %in% cell.types.ref.fail)]
  seurat.ref.obj   = SubsetData(object=seurat.ref.obj,   ident.use=cell.types.ref,   do.clean=T, subset.raw=T)
}
}

# merge the seurat objects
if(seurat.ref.addr == seurat.query.addr){
  seurat.obj = SubsetData(object=seurat.obj, ident.use=unique(c(cell.types.query, cell.types.ref)), do.clean=T, subset.raw=T)
}else{
  seurat.query.obj@meta.data$comparison_labels = eval(parse(text = sprintf("seurat.query.obj@meta.data$%s", set.ident.query)))
  seurat.query.obj@meta.data$comparison_labels = paste("Q::", as.vector(seurat.query.obj@meta.data$comparison_labels), sep = "")
  seurat.ref.obj@meta.data$comparison_labels   = eval(parse(text = sprintf("seurat.ref.obj@meta.data$%s", set.ident.ref)))
  seurat.ref.obj@meta.data$comparison_labels   = paste("R::", as.vector(seurat.ref.obj@meta.data$comparison_labels), sep = "")
  
  seurat.obj = MergeSeurat(object1=seurat.ref.obj, object2=seurat.query.obj, min.cells=0, min.genes=0, do.normalize=T, add.cell.id1="R", add.cell.id2="Q")
  seurat.obj = SetAllIdent(object=seurat.obj, id="comparison_labels")
}

# compute gene x cell type mean gene expression matrix for the query and reference data sets
query.expression.data = gene_mean_expression(seurat.query.obj)
ref.expression.data   = gene_mean_expression(seurat.ref.obj)

query.expression.data = query.expression.data[cell.types.query, ]
ref.expression.data   = ref.expression.data[cell.types.ref,     ]

# process and merge the expression matrices
if(seurat.ref.addr == seurat.query.addr){
  query.labels = cell.types.query
  ref.labels   = cell.types.ref
}else{
  query.labels = paste("Q::", rownames(query.expression.data), sep = "")
  ref.labels   = paste("R::", rownames(ref.expression.data),   sep = "")
}

rownames(query.expression.data) = query.labels
rownames(ref.expression.data)   = ref.labels

# must compute variable genes in both seurat objects, otherwise the might be no shared variable genes
print("Computing variable genes ... ")
seurat.query.obj = NormalizeData(object = seurat.query.obj, normalization.method = "LogNormalize", scale.factor = 10000)
print("Computing variable genes ... ")
seurat.query.obj = FindVariableGenes(object = seurat.query.obj, mean.function = ExpMean, 
                                     dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                     x.high.cutoff = 3, y.cutoff = .625)
seurat.ref.obj = NormalizeData(object = seurat.ref.obj, normalization.method = "LogNormalize", scale.factor = 10000)
print("Computing variable genes ... ")
seurat.ref.obj = FindVariableGenes(object = seurat.ref.obj, mean.function = ExpMean, 
                                   dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                   x.high.cutoff = 3, y.cutoff = .625)
print("Finished computing the variable genes.")

shared.genes = c(seurat.query.obj@var.genes, seurat.ref.obj@var.genes)
shared.genes = intersect(shared.genes, colnames(query.expression.data))
shared.genes = intersect(shared.genes, colnames(ref.expression.data))

query.expression.data = query.expression.data[, shared.genes]
ref.expression.data   = ref.expression.data[, shared.genes]

expression.data = rbind(query.expression.data, ref.expression.data)

# compute correlation distances
cor.dist = cor(t(expression.data), method="spearman")
cor.dist = cor.dist[query.labels, ref.labels]

p.values = c()
for(i in seq_along(1:nrow(expression.data))){
  for(j in seq_along(1:nrow(expression.data))){
    p.values = c(p.values, cor.test(x = as.numeric(expression.data[i, ]), y = as.numeric(expression.data[j, ]), method="spearman")$p.value)
  }
}
p.values = matrix(p.values, nrow=nrow(expression.data))
rownames(p.values) = rownames(expression.data)
colnames(p.values) = rownames(expression.data)
p.values = p.values[query.labels, ref.labels]
write.csv(p.values, file.path(output_folder, "correlation_p_values.csv"))

# make correlation plot
cor.plot = make_heatmap(table.data=cor.dist, measure.name="Correlation")
cor.plot.fname = file.path(output_folder, "correlation_matrix.pdf")
pdf(cor.plot.fname, width=dims.plot[1], height = dims.plot[2])
print(cor.plot)
dev.off()

####### compute AGA scores ######
print("Computing AGA scores")
print("... but first need to under-sample cell populations that are more than 500 in numbers in order to make a smaller seurat object.")
seurat.obj = DownSample(seurat.obj, limit = 2000)
print("Done sampling, back to AGA score computation ...")

AGA_folder = file.path(output_folder, "AGA_folder")
dir.create(AGA_folder)
# write cell labels to disk
write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)

# save raw data to disk
raw_data = seurat.obj@raw.data
raw_data = raw_data[rownames(seurat.obj@data), colnames(seurat.obj@data)]

writeMM(raw_data, file.path(output_folder_material, "raw_data.mtx"))
# save gene names
gene_names = rownames(raw_data)
write.csv(data.frame(Genes = gene_names), file.path(output_folder_material, "genenames.csv"))
# save cell names
cell_names = colnames(raw_data)
write.csv(data.frame(Cells = cell_names), file.path(output_folder_material, "cellnames.csv"))
# write cell labels to disk
write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)

command = file.path(tool_addr, "AGA/AGA_from_Seurat.py")
command = paste(paste(python.addr, command, sep = " "), output_folder, sep = " ")
command = paste(command, "AGA_scores", sep =" ")
system(command, wait = T)

# read the AGA output
connectivities = read.csv(file.path(AGA_folder, "connectivities.csv"), row.names = 1)
colnames(connectivities) = rownames(connectivities)
connectivities = connectivities[query.labels, ref.labels] 
print("printing connectivities from line 394")
print(connectivities)

# make AGA plot
AGA.plot = make_heatmap(table.data=connectivities, measure.name="AGA score", "cadetblue1", "#928c91", "brown4")
AGA.plot.fname = file.path(output_folder, "AGA_scores_matrix.pdf")
pdf(AGA.plot.fname, width = dims.plot[3], height = dims.plot[4])
print(AGA.plot)
dev.off()

# delete unused folders
unlink(output_folder_material,              recursive=T, force=T)
unlink(AGA_folder,                          recursive=T, force=T)
#unlink(file.path(output_folder, "figures"), recursive=T, force=T)

# clear up some memory
rm(seurat.query.obj, seurat.ref.obj, query.expression.data, ref.expression.data, expression.data, cor.dist, raw_data)

###### compute DEGs #######
if(compute.DEGs){
  print("Computing DEGs ...")
  deg.span = 50
  padW = 2
  padH = 2
  bitW = .6
  bitH = 1.1
  plot.W = 120 * (padW + 2 * bitW * deg.span)
  plot.H = 130 * (padH + bitH * length(ref.labels))
  deg.folder = file.path(output_folder, "DEGs")
  dir.create(deg.folder)
  
  print("Plotting DEGs")
  #bplapply(query.labels, plot_DEGs, BPPARAM=MulticoreParam(3))
  print("printing query labels")
  print(query.labels)
  for (i in 1:length(query.labels)){
    plot_DEGs(query.labels[i])
  }
  file.remove("Rplots.pdf")

  print("removed rplots")
  
  # merge all DEG files
  DEG_folder    = file.path(output_folder, "DEGs")
  print("made deg folder")
  DEGs.fnames   = list.files(DEG_folder)
  DEG.csv.file.paths = DEGs.fnames[unlist(lapply(DEGs.fnames, function(fname){substr(x=fname, start=nchar(fname) - 3, stop=nchar(fname)) == '.csv'}))]
  DEG.csv.file.paths = file.path(DEG_folder, DEG.csv.file.paths)
  print("made DEG.csv.file.paths")
  DEG.csv.files = lapply(DEG.csv.file.paths, read.csv, row.names = 1)
  DEG.csv.files = Reduce(f=rbind, x=DEG.csv.files)
  print("DEG.csv.files")
  write.csv(DEG.csv.files, file.path(DEG_folder, "DEGs.csv"))
  file.remove(DEG.csv.file.paths)
}

print("Ended beautifully ... ")

