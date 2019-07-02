args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg     = args[1]
plot.by.arg         = args[2]
type.to.colours.arg = args[3]
runDiffusionMap.arg = args[4]
runAGA.arg          = args[5]
"

plotW = 8
plotH = 8

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
output_folder = paste("11_plot_dr", seurat.addr, sep = "_")
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
library(RColorBrewer)
library(plyr)
library(dplyr)
library(destiny)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

# save raw data to disk
raw_data = seurat.obj@raw.data
raw_data = raw_data[rownames(seurat.obj@data), colnames(seurat.obj@data)]
  
# save gene names
gene_names = rownames(raw_data)

# save cell names
cell_names = colnames(raw_data)
write.csv(data.frame(Cells = cell_names), file.path(output_folder_material, "cellnames.csv"))
# write cell labels to disk
write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)

# run the diffusion map
if(runDiffusionMap){
  print("Writing .mtx file for diffusion map ...")
  writeMM(raw_data, file.path(output_folder_material, "raw_data.mtx"))
  
  print("Running the diffusion map ... ")
  command = sprintf("%s ./dm.py %s", python.addr, output_folder_material)
  system(command, wait = T)
  
  # load dm from disk
  dm = read.csv(file.path(output_folder_material, "dm.csv"), row.names = 1, header = F)
  #dm = DiffusionMap(seurat.obj@dr$pca@cell.embeddings, k = 100, density_norm=F, n_eigs = 20)
  #dm = data.frame(DC1 = dm$DC1, DC2 = dm$DC2, DC2 = dm$DC3)
}

print("Computing FDG limits ...")
fdg.x = seurat.obj@dr$fdg@cell.embeddings[, 1]
fdg.y = seurat.obj@dr$fdg@cell.embeddings[, 2]
fdg.limits = 1.15 * c(quantile(fdg.x, c(.01)), quantile(fdg.x, c(.99)), quantile(fdg.y, c(.01)), quantile(fdg.y, c(.99)))

print("Making the plots ...")
for (index in 1:length(plot.by)){
  caty = plot.by[index]
  seurat.obj = SetAllIdent(object=seurat.obj, id=caty)
  if (!is.na(type.to.colours[index])){
    type.to.colour = read.csv(file.path("../../resources", type.to.colours[index]))
    filter.key   =  type.to.colour$CellTypes %in% as.vector(unique(seurat.obj@ident))
    cell.labels  = as.vector(type.to.colour$CellTypes[filter.key])
    cell.colours = as.vector(type.to.colour$Colours[filter.key])
  }else{
    cell.labels  = sort(as.vector(unique(seurat.obj@ident)))
    cell.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
  }
  caty = gsub(pattern="\\.", replacement="_", caty)
  # file paths for annotated graphs
  tsne.file.name = file.path(output_folder, paste("tsne", paste(caty, "pdf", sep = "."), sep = "_"))
  umap.file.name = file.path(output_folder, paste("umap", paste(caty, "pdf", sep = "."), sep = "_"))
  fdg.file.name = file.path(output_folder, paste("fdg", paste(caty, "pdf", sep = "."), sep = "_"))
  legend.file.name = file.path(output_folder, paste("legend", paste(caty, "pdf", sep = "."), sep = "_"))
  AGA.file.name = file.path(output_folder, paste("AGA", paste(caty, "pdf", sep = "."), sep = "_"))
  
  # file paths for unannotated plots
  tsne.file.name_unlabeled = file.path(output_folder, paste("tsne_unlabeled", paste(caty, "pdf", sep = "."), sep = "_"))
  umap.file.name_unlabeled = file.path(output_folder, paste("umap_unlabeled", paste(caty, "pdf", sep = "."), sep = "_"))
  fdg.file.name_unlabeled = file.path(output_folder, paste("fdg_unlabeled", paste(caty, "pdf", sep = "."), sep = "_"))
  
  # preparing data frame
  df = data.frame(CellType = as.vector(seurat.obj@ident),
                   tSNEx = seurat.obj@dr$tsne@cell.embeddings[, 1],
                   tSNEy = seurat.obj@dr$tsne@cell.embeddings[, 2],
                   UMAPx = seurat.obj@dr$umap@cell.embeddings[, 1],
                   UMAPy = seurat.obj@dr$umap@cell.embeddings[, 2],
                   FDGx  = seurat.obj@dr$fdg@cell.embeddings[, 1],
                   FDGy  = seurat.obj@dr$fdg@cell.embeddings[, 2])
  
  interactive_plot_df = data.frame(X = seurat.obj@dr$tsne@cell.embeddings[, 1],
                                    Y = seurat.obj@dr$tsne@cell.embeddings[, 2])
  interactive_plot_df$Labels = factor(seurat.obj@ident, levels = cell.labels)
  interactive_plot_df$Colours = mapvalues(x = interactive_plot_df$Labels, from = cell.labels, to = cell.colours)
  
  # make interartive tsne 
  interactive_tsne_filename = file.path(output_folder, paste(paste("Interactive_tSNE", caty, sep = "_"), "html", sep = "."))
  make_2D_interactive_page(data_frame_2D=interactive_plot_df, tool_addr=tool_addr, python.addr=python.addr, save.to=interactive_tsne_filename)
  
  # make interactive UMAP
  interactive_plot_df$X = seurat.obj@dr$umap@cell.embeddings[, 1]
  interactive_plot_df$Y = seurat.obj@dr$umap@cell.embeddings[, 2]
  interactive_umap_filename = file.path(output_folder, paste(paste("Interactive_UMAP", caty, sep = "_"), "html", sep = "."))
  make_2D_interactive_page(data_frame_2D=interactive_plot_df, tool_addr=tool_addr, python.addr=python.addr, save.to=interactive_umap_filename)
  
  # make interactive FDG
  interactive_plot_df$X = seurat.obj@dr$fdg@cell.embeddings[, 1]
  interactive_plot_df$Y = seurat.obj@dr$fdg@cell.embeddings[, 2]
  interactive_fdg_filename = file.path(output_folder, paste(paste("Interactive_FDG", caty, sep = "_"), "html", sep = "."))
  make_2D_interactive_page(data_frame_2D=interactive_plot_df, tool_addr=tool_addr, python.addr=python.addr, save.to=interactive_fdg_filename)
  
  n.cols = min(2, length(cell.labels))
  n.rows = ceiling(length(cell.labels) / n.cols)
  
  # making the plots
  print("making the plots")
  
  # annotated plots
  plot.tsne   = dr.plot.numerical(point.labels=df$CellType, dr1=df$tSNEx, dr2=df$tSNEy, dr1.name="tSNE-x", dr2.name="tSNE-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, use.labels=cell.labels)
  plot.umap   = dr.plot.numerical(point.labels=df$CellType, dr1=df$UMAPx, dr2=df$UMAPy, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, use.labels=cell.labels)
  plot.fdg    = dr.plot.numerical(point.labels=df$CellType, dr1=df$FDGx, dr2=df$FDGy, dr1.name="FDG-x", dr2.name="FDG-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, limits=fdg.limits, use.labels=cell.labels)
  #plot.legend = plot.indexed.legend(label.vector=cell.labels, color.vector=cell.colours, ncols=n.cols, left.limit=.2, symbol.size=10, text.size=6, padH=.6, padV=.6)
  
  # unannotated plots
  plot.tsne_unlabeled   = dr.plot.numerical(point.labels=df$CellType, dr1=df$tSNEx, dr2=df$tSNEy, dr1.name="tSNE-x", dr2.name="tSNE-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, use.labels=cell.labels, annotate.plot = F)
  plot.umap_unlabeled   = dr.plot.numerical(point.labels=df$CellType, dr1=df$UMAPx, dr2=df$UMAPy, dr1.name="UMAP-x", dr2.name="UMAP-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, use.labels=cell.labels, annotate.plot = F)
  plot.fdg_unlabeled    = dr.plot.numerical(point.labels=df$CellType, dr1=df$FDGx, dr2=df$FDGy, dr1.name="FDG-x", dr2.name="FDG-y", no.legend=T, plt.lb.sz=7, txt.lb.size=4, use.cols=cell.colours, limits=fdg.limits, use.labels=cell.labels, annotate.plot = F)
  
  # print the annotated plots
  pdf(tsne.file.name, width = plotW, height = plotH)
  print(plot.tsne)
  dev.off()
  
  pdf(umap.file.name, width = plotW, height = plotH)
  print(plot.umap)
  dev.off()
  
  pdf(fdg.file.name, width = plotW, height = plotH)
  print(plot.fdg)
  dev.off()
  
  #pdf(legend.file.name, width = 1.5 + .15 * n.cols * max(unlist(lapply(cell.labels, nchar))), height = .5 + n.rows * .35)
  #print(plot.legend)
  #dev.off()
  
  # print the unannotated plots
  pdf(tsne.file.name_unlabeled, width = plotW, height = plotH)
  print(plot.tsne_unlabeled)
  dev.off()
  
  pdf(umap.file.name_unlabeled, width = plotW, height = plotH)
  print(plot.umap_unlabeled)
  dev.off()
  
  pdf(fdg.file.name_unlabeled, width = plotW, height = plotH)
  print(plot.fdg_unlabeled)
  dev.off()
  
  # run diffusion map
  if(runDiffusionMap){
    df = as.data.frame(dm[, 1:3])
    df$Labels = factor(seurat.obj@ident, levels = cell.labels)
    df$Colours = mapvalues(x = df$Labels, from = cell.labels, to = cell.colours)
    
    dm.file.name = file.path(output_folder_material, paste(paste("dm_data", caty, sep="_"), "csv", sep = "."))
    write.csv(df, dm.file.name, row.names = F)
    
    dm.file.name = file.path(output_folder, paste(paste("DiffusionMap_3D", caty, sep = "_"), "html", sep = "."))
    make_3D_interactive_page(data_frame_3D=df, tool_addr=tool_addr, python.addr=python.addr, save.to=dm.file.name)
  }
  
  if(runAGA){
    if(runDiffusionMap==F){
      print("Writing .mtx file for AGA map ...")
      writeMM(raw_data, file.path(output_folder_material, "raw_data.mtx"))
    }
    print("running AGA ...")
    AGA_folder = file.path(output_folder, "AGA_folder")
    dir.create(AGA_folder)
    # write cell labels to disk
    write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)
    
    # running AGA
    command =file.path(tool_addr, "AGA/AGA_from_Seurat.py")
    command = paste(paste(python.addr, command, sep = " "), output_folder, sep = " ")
    command = paste(command, caty, sep =" ")
    system(command, wait = T)
    
    # read the AGA output
    coordinates = read.csv(file.path(AGA_folder, "coordinates.csv"), row.names = 1)
    connectivities = read.csv(file.path(AGA_folder, "connectivities.csv"), row.names = 1)
    
    # plot AGA
    coordinates = coordinates[cell.labels, ]
    coordinates$Colours = cell.colours
    label.order = match(cell.labels, rownames(connectivities))
    connectivities = connectivities[label.order, label.order]
    
    plot.obj = ggplot(data=coordinates, aes(x = X, y = Y))
    plot.obj = plot.obj + theme_void() + theme(legend.position="none")
    xi = c(); xf = c(); yi = c(); yf = c(); vs = c();
    for(i in 1:dim(connectivities)[1]){
      for(j in i:dim(connectivities)[2]){
        v = connectivities[i, j]
        if(v > 0){
          xi = c(xi, coordinates$X[i])
          xf = c(xf, coordinates$X[j])
          yi = c(yi, coordinates$Y[i])
          yf = c(yf, coordinates$Y[j])
          vs = c(vs, v)
        }
      }
    }
    lineDF = data.frame(Xi = xi, Yi = yi, Xf = xf, Yf = yf, Vs = vs)
    plot.obj = plot.obj + geom_segment(data = lineDF, aes(x = Xi, y = Yi, xend = Xf, yend = Yf), color = "black", size = 3 * lineDF$Vs)
    plot.obj = plot.obj + geom_point(size = 2 * log(coordinates$Size), color = coordinates$Colours)
    plot.obj = plot.obj + annotate("text", x=coordinates$X, y=coordinates$Y, label = 1:dim(coordinates)[1])
    pdf(AGA.file.name, width = 10, height = 10)
    print(plot.obj)
    dev.off()
    
    ######## now make the interactive AGA app 
    #########################################
    print("Making the AGA app ... ")
    # save colours
    colours.df = data.frame(CellTypes = cell.labels, Colours = cell.colours)
    write.csv(colours.df, file.path(AGA_folder, "colours.csv"), row.names = F)
    
    # run python to built the AGA app
    command = sprintf("%s make_AGA_app.py %s %s", python.addr, output_folder, caty)
    system(command, wait = T)
  }
}

# cleaning garbage folders
unlink(output_folder_material, recursive=T, force=T)
if (runAGA){
  unlink(AGA_folder,                          recursive=T, force=T)
  unlink(file.path(output_folder, 'figures'), recursive=T, force=T) 
}

print("Ended beautifully ... ")
