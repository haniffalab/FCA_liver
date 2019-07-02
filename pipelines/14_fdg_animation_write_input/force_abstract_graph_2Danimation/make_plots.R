setwd("~/Documents/MyTools/force_abstract_graph_2Danimation/")
buffers.addrs <- list.files("./input/buffers/", full.names=T)
data.colours <- as.vector(read.csv("./input/label_colours.csv")$LabelCols)

################################################################################################################
################################################################################################################
################################################################################################################
library(RColorBrewer)
library(dplyr)
library(plyr)
library(Seurat)

#c.unique  <-as.vector( unique(data.colours))
#c.colours <- sample(colorRampPalette(brewer.pal(12, "Paired"))(length(c.unique)))
#data.colours <- factor(plyr::mapvalues(x=data.colours, from=c.unique, to = c.colours), levels = c.colours)

################################################################################################################
################################################################################################################
################################################################################################################

for(k in 1:length(buffers.addrs)){
  buffer.addr <- buffers.addrs[k]
  print(sprintf("Plotting frame %d", k))
  buffer.data <- read.csv(buffer.addr, header = F)
  buffer.data <- cbind(buffer.data, data.colours)
  colnames(buffer.data) <- c("FDGX", "FDGY", "Colours")
  limitX <- quantile(buffer.data$FDGX, c(.01, .99)) + c(-15000, 15000)
  limitY <- 1.1 * quantile(buffer.data$FDGY, c(.01, .99)) + c(-15000, 15000)
  plot.obj <- ggplot(data=buffer.data, aes(x = FDGX, y = FDGY))
  plot.obj <- plot.obj + geom_point(show.legend=F, size = 1.5, color = as.vector(buffer.data$Colours))
  plot.obj <- plot.obj + scale_color_manual(values=as.vector(buffer.data$Colours))
  plot.obj <- plot.obj + theme(plot.background = element_rect(fill = "black"))
  plot.obj <- plot.obj + scale_x_continuous(limits = limitX, expand = c(0, 0))
  plot.obj <- plot.obj + scale_y_continuous(limits = limitY, expand = c(0, 0))
  plot.obj <- plot.obj + theme(axis.title = element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank())
  fname <- file.path("./input/frames", sub(pattern=".csv", replacement=".png", x=basename(buffer.addr)))
  png(fname, width = 2000, height = 2000)
  print(plot.obj)
  dev.off()
}
