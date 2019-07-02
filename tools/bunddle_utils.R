tool_addr = "../../tools"
python.addr = "python3.6"

# a function to compute force-directed graph; return coordinates as a data frame
# arguments:
# pca.df: pca data as a data frame
# snn   : a nearest-neighbor graph as a sparse data matrix
runFDG = function(pca.df, snn, iterations = 600, tool_addr, python.addr){
  current.wd = getwd()
  setwd(file.path(tool_addr, "force_abstract_graph_2D"))
  # generate unique name for pca data file
  pca.data.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  pca.data.fname = paste(pca.data.fname, ".csv", sep = "")
  # generate unique name for snn file
  snn.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  snn.fname = paste(snn.fname, ".smm", sep = "")
  # generate unique name for fdg coordinates
  fdg.coordinates.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  fdg.coordinates.fname = paste(fdg.coordinates.fname, ".csv", sep = "")
  write.csv(pca.df, pca.data.fname)
  writeMM(obj=snn, file=snn.fname)
  command = gsub(pattern="ITER", replacement=as.character(iterations), paste(python.addr, "./make_fdg.py ITER", sep = " "))
  command = paste(command, paste(c(pca.data.fname, snn.fname, fdg.coordinates.fname), collapse = " "), sep = " ")
  system(command, wait = T)
  fdg_coordinates = read.csv(fdg.coordinates.fname, header = FALSE)
  colnames(fdg_coordinates) = c("X", "Y")
  rownames(fdg_coordinates) = rownames(pca.df)
  file.remove(c(pca.data.fname, snn.fname, fdg.coordinates.fname))
  setwd(current.wd)
  return(fdg_coordinates)
}

# a function to perform UMAP on a seurat object using PCA coordinates
RunUMAP = function(pca.df, tool_addr, python.addr){
  current.wd = getwd()
  setwd(file.path(tool_addr, "umap"))
  print("writting pca data to disk...")
  # generate unique name for pca data file
  pca.data.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  pca.data.fname = paste(pca.data.fname, ".csv", sep = "")
  # generate unique name for umap coordinates data file
  umap.coordinates.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  umap.coordinates.fname = paste(umap.coordinates.fname, ".csv", sep = "")
  write.csv(pca.df, pca.data.fname)
  print("perfoming UMAP")
  command = paste(python.addr, "umap_compute.py", sep = " ")
  command = paste(command, pca.data.fname, sep = ' ')
  command = paste(command, umap.coordinates.fname, ' ')
  system(command, wait = T)
  print("Reading results...")
  umap.coordinates = read.csv(umap.coordinates.fname, stringsAsFactors = F)
  file.remove(c(pca.data.fname, umap.coordinates.fname))
  setwd(current.wd)
  umap.coordinates = umap.coordinates[, c("UMAPx", "UMAPy")]
  return(umap.coordinates)
}

# a function that load a SVM models and make predictions like cell types or doublet.singlets
Apply_Classifier_On_Seurat_Object = function(seurat.obj, classifier.fname, tool_addr, python.addr){
  current.wd = getwd()
  setwd(file.path(tool_addr, 'predict_by_classifier'))
  predictor.addr = file.path("../../resources", classifier.fname, sep = "")
  print(predictor.addr)
  if(!dir.exists(predictor.addr)){
    print("classifier does not exists")
    available.classifiers = list.dirs("../../resources", full.names=F)
    available.classifiers = available.classifiers[grepl("classifier_", available.classifiers)]
    print(paste("Available doublets identifiers: ", paste(available.classifiers, collapse = ", "), sep = ""))
    setwd(current.wd)
    return(NULL)
  }
  tryCatch({
    OK = FALSE
    print("reading feature genes ...")
    feature.genes = readRDS(file.path(predictor.addr, "feature_genes.RDS"))
    features.present = feature.genes[feature.genes %in% rownames(seurat.obj@data)]
    features.not.present = feature.genes[!(feature.genes %in% rownames(seurat.obj@data))]
    expr.data = as.data.frame(t(as.matrix(seurat.obj@data[features.present, ])))
    if (length(features.not.present) > 0){
      zeros.data = as.data.frame(t(as.matrix(seurat.obj@data[1:length(features.not.present),])))
      zeros.data[] = 0
      colnames(zeros.data) = features.not.present
      expr.data=cbind(expr.data, zeros.data)
      expr.data = expr.data[, feature.genes]
    }
    print("Writting data to disk ... ")
    data.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
    data.fname = paste(data.fname, ".csv", sep = '')
    write.csv(expr.data, data.fname)
    
    print(sprintf("Dims are: %s", dim(expr.data)))
    print("Copying pickle file ... ")
    model.fname = paste(sample(LETTERS, 20, T), collapse = '')
    model.fname = paste(model.fname, '.pickle', sep = "")
    file.copy(from=file.path(predictor.addr, "model.pickle"), to=model.fname)
    
    pca.fname = paste(sample(LETTERS, 20, T), collapse = '')
    pca.fname = paste(pca.fname, '.pickle', sep = "")
    file.copy(from=file.path(predictor.addr, "pca.pickle"), to=pca.fname)
    
    print("Running classifier in Python ... ")
    predictions.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
    predictions.fname = paste(predictions.fname,  ".csv", sep = "")
    command = paste(python.addr, "predict.py", sep = " ")
    command = paste(command, paste(c(model.fname, data.fname, predictions.fname, pca.fname), collapse = " "), sep = " ")
    system(command, wait = T)
    predictions = read.csv(predictions.fname)
    OK = TRUE
  },
  warning = function(warning_conditions){print("")},
  error = function(error_condition){
    print("Errors occured. Cleaning up and then returning NULL ...")
    file.remove(c(model.fname, data.fname, predictions.fname, pca.fname))
    setwd(current.wd)
  })
  if(OK){
    if(length(features.not.present) > 0){
      print('Everything went well except the fact that some of the features genes were not present in the data. These are:')
      print(paste(features.not.present, collapse = ", "))
    }else{
      print("Everything went smooth. Cleaning up and returning the predictions ... ")
    }
    file.remove(c(model.fname, data.fname, predictions.fname, pca.fname))
    setwd(current.wd)
    return(predictions$X0)
  }
  return(NULL)
}

# a function that takes 3D coordinates (e.g. from diffusion map) and generates an
# interactive html page
make_3D_interactive_page = function(data_frame_3D, tool_addr, python.addr, save.to){
  data.frame.fname = paste(sample(LETTERS, 20, T), collapse = '')
  data.frame.fname = paste(data.frame.fname, ".csv", sep = "")
  save.to = file.path(getwd(), save.to)
  command = gsub(pattern="save.to", replacement=save.to, x=paste(python.addr, './html_WebGL_3D_viewer.py save.to', sep = " "))
  command = paste(command, data.frame.fname, sep = " ")
  current.wd = getwd()
  setwd(file.path(tool_addr, "interactive_3D_viewer"))
  write.csv(data_frame_3D, data.frame.fname, row.names = F)
  system(command, wait=T)
  file.remove(data.frame.fname)
  setwd(current.wd)
}

# a function that takes 2D coordinates and generates an interactive html page
make_2D_interactive_page = function(data_frame_2D, tool_addr, python.addr, save.to="./"){
  data.frame.fname = paste(sample(LETTERS, 20, T), collapse = '')
  data.frame.fname = paste(data.frame.fname, ".csv", sep = "")
  save.to = file.path(getwd(), save.to)
  command = gsub(pattern="save.to", replacement=save.to, x=paste(python.addr, './html_WebGL_2D_viewer.py save.to', sep = " "))
  command = paste(command, data.frame.fname, sep = " ")
  current.wd = getwd()
  setwd(file.path(tool_addr, "interactive_2D_viewer"))
  write.csv(data_frame_2D, data.frame.fname, row.names = F)
  system(command, wait=T)
  file.remove(data.frame.fname)
  setwd(current.wd)
}

# a function that create interactic html pages to explore gene expression in data parsed by different categories
create_gene_expression_viewer_apps = function(seurat.obj, dim.type = 'umap', save.to, tool_addr, python.addr, categories.colours=NA){
  categories = c("cell.labels", "fetal.ids", "sort.ids", "lanes", "stages", "gender", "doublets")
  if(is.na(categories.colours)){
    categories.colours = rep(NA, length(categories))
  }
  categories.data = as.data.frame(seurat.obj@meta.data[names(seurat.obj@ident), categories])
  for(j in 1:length(categories)){
    category = categories[j]
    category.colour.scheme = categories.colours[j]
    if (!is.na(category.colour.scheme)){
      category.colour.scheme = read.csv(category.colour.scheme)
      category.colour.scheme = mapvalues(x=categories.data[, category], from=as.vector(unique(category.colour.scheme$CellTypes)), to=as.vector(unique(category.colour.scheme$Colours)))
    }else{
      category.colour.scheme = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(as.vector(unique(categories.data[, category])))))
      category.colour.scheme = mapvalues(x=categories.data[, category], from=as.vector(unique(categories.data[, category])), to=category.colour.scheme)
    }
    category = paste(category, "colours", sep = "_")
    categories.data[, category] = category.colour.scheme
  }
  eval(parse(text = sprintf("dim.data = seurat.obj@dr$%s@cell.embeddings[names(seurat.obj@ident), 1:2]", dim.type)))
  n.categories = length(categories)
  
  genes.by.file = round( 20 * 26149846 / ncol(seurat.obj@data))
  approx.no.of.files = round(nrow(seurat.obj@data) / genes.by.file)
  gene.names = sort(rownames(seurat.obj@data))
  gene.names = sort(seurat.obj@var.genes)
  gene.splits = split(gene.names, sort(1:length(gene.names) %% approx.no.of.files))
  curdir = getwd()
  unlink(x=save.to, recursive=T)
  dir.create(save.to)
  setwd(file.path(tool_addr, 'gene_expression_viewer_apps'))
  folder_name = paste(sample(LETTERS, 20, T), collapse = "")
  dir.create(folder_name)
  for (l in 1:length(gene.splits)){
    gene.split = unlist(gene.splits[[l]])
    expression.data = as.data.frame(as.matrix(t(seurat.obj@data[gene.split, names(seurat.obj@ident)])))
    expression.data = cbind(dim.data, categories.data, expression.data)
    first_gene = gene.split[1]
    last_gene = gene.split[length(gene.split)]
    print(sprintf("Creating app for: %s", paste(c(first_gene, "to", last_gene), collapse = "_")))
    expression.data.fname = paste(c(first_gene, "to", last_gene), collapse = "_")
    expression.data.fname = paste(expression.data.fname, ".csv", sep = "")
    expression.data.fname = file.path(folder_name, expression.data.fname)
    save_to = paste(c(first_gene, "to", last_gene), collapse = "_")
    save_to = paste(save_to, ".html", sep = "")
    save_to = file.path(file.path(curdir, save.to), save_to)
    command = sprintf("%s gene_expression_viewer_apps.py %s %s %s", python.addr, save_to, expression.data.fname, n.categories)
    write.csv(expression.data, expression.data.fname, row.names = F)
    system(command, wait = T)
  }
  unlink(x=folder_name, recursive=T, force=T)
  setwd(curdir)
}

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
dr.plot = function(point.labels, dr1, dr2, dr1.name, dr2.name, no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2, use.cols = NULL, use.labels = NULL, limits = NULL, annotate.plot = T, index.map = NA){
  if(!is.na(overlay.data)){
  df.dr = data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2, overlay.data=factor(df$overlay.data))
  }
  else{
  df.dr = data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2)
  }
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
  if(is.na(index.map)){
    index.map = 1:length(p.labels)
  }
  if(!is.na(overlay.data)){
  plot.obj = ggplot(data = df.dr, aes(x = DR1, y = DR2, color = Cell.Labels, shape = factor(overlay.data, levels=overlay.data.ordered)))
  print("levels=overlay.data.ordered")
  print(overlay.data.ordered)
  }
  else{
  plot.obj = ggplot(data = df.dr, aes(x = DR1, y = DR2, color = Cell.Labels))
  }
  # this line should give different sizes for different values in overlay.data metadata column
  if(!is.na(overlay.data)){
  plot.obj = plot.obj + geom_point(size=pt.size)
  plot.obj = plot.obj + geom_point(data = subset(df.dr, overlay.data == overlay.data.ordered[2]))
  print("overlay.data.ordered[2]")
  print(overlay.data.ordered[2])
  }
  else{
  plot.obj = plot.obj + geom_point(size=pt.size)
  }
  plot.obj = plot.obj + scale_color_manual(values=plt.colours)
  if(annotate.plot){
    if(!is.na(overlay.data)){
    plot.obj = plot.obj + geom_point(data=p.labels.medians,aes(x = DR1, y = DR2), colour = "gray", size = plt.lb.sz, fill = plt.colours, alpha = .5, pch = 21, shape=factor(df.dr$overlay.data, levels=overlay.data.ordered))
    }
    else{
    plot.obj = plot.obj + geom_point(data=p.labels.medians,aes(x = DR1, y = DR2), colour = "gray", size = plt.lb.sz, fill = plt.colours, alpha = .5, pch = 21)
    }
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
                    

# plotting function for dimensionaly-reduced data to label population by a round indexed label
dr.plot.numerical = function(point.labels, dr1, dr2, dr1.name, dr2.name, no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random_state = 2, use.cols = NULL, use.labels = NULL, limits = NULL, annotate.plot = T, index.map = NA){
  df.dr = data.frame("Cell Labels" = point.labels, DR1 = dr1, DR2 = dr2)
  if(is.null(use.labels)){
    p.labels = sort(unique(as.vector(point.labels)))
  }
  else{
    p.labels = use.labels
  }
  df.dr$Cell.Labels = factor(df.dr$Cell.Labels, levels=p.labels)
  p.labels.medians = aggregate(df.dr[, 2:3], list(df.dr$Cell.Labels), median)
  #df.dr$Cell.Labels = mapvalues(x = df.dr$Cell.Labels, from = p.labels, to = paste(1:length(p.labels), p.labels, sep = " "))
  if(is.null(use.cols)){
    set.seed(random_state)
    plt.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(p.labels)))
  }else{
    plt.colours = use.cols
  }
  if(is.na(index.map)){
    #index.map = 1:length(p.labels)
    index.map = p.labels
  }
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
