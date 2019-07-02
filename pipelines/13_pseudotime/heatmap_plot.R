# Prepare a smaller pseudotime heatmap, using the following genes:
selected.gene.list <- scan("selected.genes.std.txt", what = character(), sep = "\n", blank.lines.skip = T, comment.char = "#") # or character vector c("")
path <- "." # path to 'ploting.material.RDS' [sic]
library("ggplot2")
###############################################################################
plottingmat <- readRDS(file.path(path, "ploting_material.RDS"))
# str(plottingmat)
# str(plottingmat$beautiful_result_norm)
# View(plottingmat$beautiful_result_norm)

subsetplotmat <- plottingmat$beautiful_result_norm[plottingmat$beautiful_result_norm$GeneNames %in% selected.gene.list, ]
subsetplotmat$GeneNames <- droplevels(subsetplotmat$GeneNames)
subsetplotmat$GeneNames <- factor(subsetplotmat$GeneNames, levels = rev(selected.gene.list)) # Orders the heatmap

# The following section is adapted from: https://github.com/haniffalab/Single-cell-RNAseq-data-analysis-bundle/blob/master/pipelines/13_pseudotime/pseudotime.R#L270 commit b86d20dc87d35820daac178a93e46badf99216ab
plot.genes <- ggplot(data = subsetplotmat, aes(x = Pseudotime, y = GeneNames))
plot.genes <- plot.genes + geom_tile(aes(fill = ExpressionValue),
  width=1.001, height=1.001)
plot.genes <- plot.genes + scale_fill_gradient2(low = "deepskyblue",
  high = "firebrick3",
  mid = "darkolivegreen3",
  midpoint = 0.5,
  name = "Minmax normalized gene expression")
plot.genes <- plot.genes + theme(legend.position = "bottom",
  legend.text = element_text(size = 25, angle = 90),
  legend.title = element_text(size = 25),
  legend.key.width = unit(2, "cm"),
  axis.text.x = element_blank(), axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 8))

plot.genes

height = 6; width = 3
pdf("dpt_heatmap.pdf", height = height, width = width)
plot.genes
dev.off()

svg("dpt_heatmap.svg", height = height, width = width)
plot.genes
dev.off()

postscript("dpt_heatmap.ps", height = height, width = width)
plot.genes
dev.off()

png("dpt_heatmap.png", height = 600, width = 300)
plot.genes
dev.off()

###############################################################################
# Alternative formats for density plots
pdt_exp <- read.csv(file.path(path, "pdt_and_expression.csv"))
#~ str(pdt_exp)
# Standard:
ggplot(data = pdt_exp, aes(x = Pseudotime, color = Labels, fill = Labels)) + geom_density(alpha = .7) # alpha for transparency
# Stacked:
ggplot(data = pdt_exp, aes(x = Pseudotime, color = Labels, fill = Labels)) + geom_density(position = "stack")
# Relative:
ggplot(data = pdt_exp, aes(x = Pseudotime, color = Labels, fill = Labels)) + geom_density(adjust = 1.5, position = "fill")
# Histogram:
ggplot(data = pdt_exp, aes(x = Pseudotime, color = Labels, fill = Labels)) + geom_histogram(binwidth = 0.01)
