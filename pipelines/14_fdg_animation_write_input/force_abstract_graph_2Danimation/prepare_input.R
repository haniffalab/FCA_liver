# import libraries
library(Seurat)
library(plyr)

seurat.obj.addr <- "../../seurat_data/liver_immune.RDS"

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

# load the seurat object
print("Loading the data ... ")
seurat.obj <- readRDS(seurat.obj.addr)
cell.type.to.colour <- read.csv("./liver_cell_type_colours.csv")

seurat.obj <- SetAllIdent(object=seurat.obj, id="cell.labels")

################################################
print("saving pca data ...")
pca.data <- seurat.obj@dr$pca@cell.embeddings
write.csv(pca.data, "./input/pca.csv")

################################################
print("Computing and saving KNN graph ...")
seurat.obj <- BuildSNN(object=seurat.obj, reduction.type="pca", dims.use=1:20, plot.SNN=F, force.recalc=T)
writeMM(obj=seurat.obj@snn, file="./input/SNN.smm")

labels        <- as.vector(seurat.obj@ident)
labels.unique <- unique(labels)
filter.key    <- cell.type.to.colour$CellTypes %in% labels.unique
cell.labels   <- cell.type.to.colour$CellTypes[filter.key]
cell.colours  <- cell.type.to.colour$Colours[filter.key]

labels.cols <- mapvalues(x=labels, from=as.vector(cell.labels), to=as.vector(cell.colours))
write.csv(data.frame(LabelCols = labels.cols), "./input/label_colours.csv")

png("./input/legend.png", width = 1000, height = 800)
legend.plt <- plot.indexed.legend(label.vector=cell.labels, color.vector=cell.colours, left.limit=3.6, text.size=10, ncols=2, symbol.size = 15)
print(legend.plt)
dev.off()

print("End")







