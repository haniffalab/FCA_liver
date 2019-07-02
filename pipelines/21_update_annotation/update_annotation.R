args = commandArgs(trailingOnly=T)

args = paste(args, collapse = "")

args = unlist(strsplit(args, ";"))




arguments.list = "

seurat.addr.arg = args[1]

make.app        = args[2]

update.file        = args[3]

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




seurat.addr = file.path("../../data", seurat.addr)




source("../../tools/bunddle_utils.R")




library(Seurat)

library(plyr)

library(dplyr)

library(reshape2)

library(RColorBrewer)

library(wordcloud)




gene_to_weighted_cell_mention = function(gene.expr){

  idx = which(as.vector(gene_to_pop$V1) %in% names(gene.expr))

  gene.expr = gene.expr[as.vector(gene_to_pop$V1)[idx]]

  pop.expr = c()

  pop.names = c()

  for (k in 1:length(idx)){

    gene.name = names(gene.expr)[k]

    gene.value = gene.expr[k]

    pop.flags = as.vector(gene_to_pop$V2)[as.vector(gene_to_pop$V1) == gene.name]

    pop.flags = unlist(strsplit(pop.flags, ", "))

    for (p in 1:length(pop.flags)){

      pop.flag = pop.flags[p]

      gene.v = 100 * gene.value / populations.weight[pop.flag]

      if (pop.flag %in% pop.names){

        pop.expr[pop.flag] = pop.expr[pop.flag] + gene.v

      }else{

        pop.names = c(pop.names, pop.flag)

        pop.expr = c(pop.expr, gene.v)

        names(pop.expr) = pop.names

      }

    }

  }

  pop.expr

}




# load data

print("loading data ... ")

seurat.obj = readRDS(seurat.addr)

print("Data loaded.")



# load updated annotation

update.template = read.csv(update.file, stringsAsFactors = F, sep = '\t')

if(dim(update.template)[2] == 1){

  update.template = read.csv(update.file, stringsAsFactors = F, sep = ',')

}




# update cell labels in seurat object

seurat.obj@meta.data$cell.labels = mapvalues(as.vector(seurat.obj@meta.data$LouvainClustering), from = update.template$Cluster, to = update.template$Identity)




print("Saving seurat object")

saveRDS(seurat.obj, seurat.addr)




if (make.app){

  print('Making the interactive app')

  marker.genes.top = read.csv("annotation_markers.csv", stringsAsFactors = F)

  # update cluster names in marker.genes.top

  marker.genes.top$cluster = mapvalues(x=as.vector(marker.genes.top$cluster), from = update.template$Cluster, to = update.template$Identity)

  

  # now make an interactive maps

  gene_sym_to_marker = marker.genes.top[, c('gene', 'cluster')]

  categories = c("LouvainClustering", "fetal.ids", "sort.ids", "lanes", "stages", "gender", "doublets", "cell.labels")

  genes = as.vector(unique(gene_sym_to_marker$gene))

  

  expression.data = as.data.frame(as.matrix(t(seurat.obj@data[genes, names(seurat.obj@ident)])))

  categories.colours = rep(NA, length(categories))

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

  dim.data = seurat.obj@dr$umap@cell.embeddings[, 1:2]

  expression.data = cbind(dim.data, categories.data, expression.data)

  write.csv(expression.data, "./expression_data.csv", row.names = F)

  

  #gene.families = as.vector(unique(unlist(strsplit(as.vector(gene_sym_to_marker$cluster), "\\|"))))

  gene_sym_to_marker$ClusterName = as.character(gene_sym_to_marker$cluster)

  #gene_sym_to_marker$ClusterName = paste('000', gene_sym_to_marker$ClusterName, sep = '')

  #gene_sym_to_marker$ClusterName = unlist(lapply(gene_sym_to_marker$ClusterName, function(cluster_name){substr(cluster_name, nchar(cluster_name) - 2, nchar(cluster_name))}))

  #gene_sym_to_marker$ClusterName = paste('Cluster', gene_sym_to_marker$ClusterName, sep = '_')

  

  gene.families = as.vector(unique(gene_sym_to_marker$ClusterName))

  

  gene.to.family = c()

  for(i in 1:length(gene.families)){

    gene.family = gene.families[i]

    gene.family = gsub(pattern="\\'", replacement="", x=gene.family)

    members = as.vector(gene_sym_to_marker$gene[grep(gene_sym_to_marker$ClusterName, pattern=gene.family, value=F)])

    inline = sprintf("gene_families['%s']=", gene.family)

    members = paste(members, "'", sep = "")

    members = paste("'", members, sep = "")

    members = paste(members, collapse = ",")

    members = paste("[", members, "]", sep = "")

    inline = paste(inline, members)

    gene.to.family = c(gene.to.family, inline)

  }

  all.genes = as.vector(unique(gene_sym_to_marker$gene))

  all.genes = paste(all.genes, "'", sep = "")

  all.genes = paste("'", all.genes, sep = "")

  all.genes = paste(all.genes, collapse = ",")

  all.genes = paste("gene_families['ALL']=[", all.genes, "]", sep = "")

  gene.to.family = c(gene.to.family, all.genes)

  gene.to.family = sort(gene.to.family)

  

  gene.families.file = file('gene_families.txt', "w")

  writeLines(gene.to.family, gene.families.file)

  close(gene.families.file)

  

  save.to = file.path(output_folder, 'interactive_markers.html')

  n_categories = length(categories)

  command = sprintf('%s html_2D_gene_expression_viewer_by_gene_family.py %s %s %s', python.addr, 

                    save.to, 'expression_data.csv', n_categories)

  system(command, wait = T)

  

  file.remove(c('./expression_data.csv', './gene_families.txt'))

  

  # make annotation clouds for each cluster

  seurat.obj = SetAllIdent(object=seurat.obj, id='cell.labels')

  

  expression.data = seurat.obj@data

  mito.genes = grep(pattern="^MT-", x=rownames(expression.data))

  expression.data = expression.data[-c(mito.genes), ]

  

  gene_to_pop = read.csv("./gene_to_pop.tsv", sep = '\t', header = F)

  populations = paste(as.vector(gene_to_pop$V2), collapse = ", ")

  populations = unlist(strsplit(populations, ", "))

  populations.table = table(populations)

  populations.weight = as.vector(populations.table)

  names(populations.weight) = names(populations.table)

  

  idents = as.vector(unique(seurat.obj@ident))

  for (i  in 1:length(idents)){

    ident = idents[i]

    print(ident)

    ident = names(seurat.obj@ident)[seurat.obj@ident == ident]

    expression.data = as.matrix(seurat.obj@data[,ident])

    expression.data = rowMeans(expression.data)

    genes = names(expression.data)

    genes = genes[!(genes %in% genes[grep(pattern='^MT-', x=genes)])]

    expression.data = expression.data[genes]

    pop.expr = gene_to_weighted_cell_mention(expression.data)

    clouder = round(100 * pop.expr)

    

    fname = sprintf('%s.pdf', idents[i])

    fname = gsub(pattern="/", replacement="-", x=fname)

    fname = file.path(output_folder, fname)

    pdf(fname, width = 10, height = 10)

    wordcloud(words=names(clouder), clouder, min.freq = 1, max.words=500, 

              random.order=FALSE, rot.per=0.0, colors=brewer.pal(8, "Dark2"),

              order.color = T)

    dev.off()

  }

}




print("Ended beautifully ... ")
