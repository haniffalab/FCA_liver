# labels have been updated, should remove the part that overwrites cell labels
# must create functions that handle the formation of FDG animation:
# - data writter
# - plotting that takes dimenssion parameters

library(plyr)
library(RColorBrewer)
library(Seurat)

seurat.addr    <- "../../data/test_yolk_sac_subset.RDS"
seurat.obj <- readRDS(seurat.addr)
cell.type.to.colour <- read.csv("../../resources/test_yolk_sac_fdg_colour_key.csv")

print("Checking for doublets:")
print(table(seurat.obj@meta.data$doublets))

# a plotting function for indexed legend; special modifications for current script
plot.indexed.legend <- function(label.vector, color.vector, ncols = 2, left.limit = 3.4, symbol.size = 8, text.size = 10){
  if (length(label.vector) != length(color.vector)){
    stop("number of labels is different from number colors\nAdvice: learn to count!")
  }
  if (length(ncol) > length(label.vector)){
    stop("You cannot have more columns than labels\nSolution: Learn to count")
  }
  indices.vector <- 1:length(label.vector)
  label.no <- length(label.vector)
  nrows <- ceiling(label.no / ncols)
  legend.frame <- data.frame(X = rep(0, label.no), Y = rep(0, label.no), CS = color.vector, Txt = label.vector)
  for (i in 1:label.no){
    col.index <- floor(i / (nrows + 1)) + 1
    row.index <- 15 - ((i - 1) %% nrows + 1)
    legend.frame[i, 1] <- (col.index - 1) * 2
    legend.frame[i, 2] <- row.index
  }
  plot.obj <- ggplot(data = legend.frame, aes(x = X, y = Y))
  plot.obj <- plot.obj + geom_point(size = symbol.size, colour = color.vector)
  plot.obj <- plot.obj + scale_x_continuous(limits = c(0, left.limit)) + theme_void()
  plot.obj <- plot.obj + annotate("text", x=legend.frame$X+.1, y = legend.frame$Y, label=legend.frame$Txt, hjust = 0, size = text.size, colour = "white")
  plot.obj <- plot.obj + theme(panel.background = element_rect(fill='black'))
  return(plot.obj)
}

# a plotting function for indexed legend
plot.indexed.legend <- function(label.vector, color.vector, ncols = 2, left.limit = 3.4, symbol.size = 8, text.size = 10, padH = 1, padV = 1, padRight = 0){
  if (length(label.vector) != length(color.vector)){
    stop("number of labels is different from number colors\nAdvice: learn to count!")
  }
  if (length(ncol) > length(label.vector)){
    stop("You cannot have more columns than labels\nSolution: Learn to count")
  }
  indices.vector <- 1:length(label.vector)
  label.no <- length(label.vector)
  nrows <- ceiling(label.no / ncols)
  legend.frame <- data.frame(X = rep(0, label.no), Y = rep(0, label.no), CS = color.vector, Txt = label.vector)
  legend.frame$X <- rep(1:ncols, each=nrows)[1:nrow(legend.frame)]
  legend.frame$Y <- rep(nrows:1, times = ncols)[1:nrow(legend.frame)]
  Xrange <- range(legend.frame$X)
  Yrange <- range(legend.frame$Y)
  plot.obj <- ggplot(data = legend.frame, aes(x = X, y = Y))
  plot.obj <- plot.obj + geom_point(size = symbol.size, colour = color.vector)
  plot.obj <- plot.obj + scale_x_continuous(limits = c(Xrange[1] - padRight, Xrange[2] + padH))
  plot.obj <- plot.obj + scale_y_continuous(limits = c(Yrange[1] - padV, Yrange[2] + padV))
  plot.obj <- plot.obj + theme_void()
  
  plot.obj <- plot.obj + annotate("text", x=legend.frame$X, y = legend.frame$Y, label = indices.vector, size = text.size)
  plot.obj <- plot.obj + annotate("text", x=legend.frame$X+.1, y = legend.frame$Y, label=legend.frame$Txt, hjust = 0, size = text.size, colour = "white")
  plot.obj <- plot.obj + theme(panel.background = element_rect(fill='black'))
  return(plot.obj)
}

pca.data <- seurat.obj@dr$pca@cell.embeddings
write.csv(pca.data, "./input/pca.csv")

seurat.obj <- BuildSNN(object=seurat.obj, reduction.type="pca", dims.use=1:20, plot.SNN=F,force.recalc=T)
writeMM(obj=seurat.obj@snn, file="./input/SNN.smm")

labels        <- as.vector(seurat.obj@meta.data$cell.labels)
labels.unique <- unique(labels)

print("printing cell.type.to.colour")
print(cell.type.to.colour)

print("!is.na(cell.type.to.colour)")
print(!is.na(cell.type.to.colour))

if(!is.na(cell.type.to.colour)){
  cell.labels  <- as.vector(cell.type.to.colour$CellTypes)
  cell.colours <- as.vector(cell.type.to.colour$Colours)
  filter.key <- cell.labels %in% labels.unique
  cell.labels <- cell.labels[filter.key]
  cell.colours <- cell.colours[filter.key]
}else{
  cell.labels <- labels.unique
  set.seed(100)
  cell.colours <- sample(colorRampPalette(brewer.pal(12, "Paired"))(length(labels.unique)))
}

print("printing cell.labels")
print(cell.labels)
print("printing cell.colours")
print(cell.colours)

labels.cols <- mapvalues(x=labels, from=cell.labels, to=cell.colours)
write.csv(data.frame(LabelCols = labels.cols), "./input/label_colours.csv")

png("./input/legend.png", width = 1000, height = 700)
legend.plt <- plot.indexed.legend(label.vector=cell.labels, color.vector=cell.colours, ncols=2, left.limit=0, symbol.size=17, text.size=10, padH=.9, padV=.6)
print(legend.plt)
dev.off()

print("ended beautifully")

