args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))
if(length(args) != 5){
  stop('This pipeline requires 5 parameters: seurat.addrs (list of file names)\nappend_tag (boolean)\ntags_to_append (list of tags)\nappend_tags_at (list of meta.data columns where to append the tags)\nsave (file name to save the data at)')
}

arguments.list = "
seurat.addrs.arg    = args[1]
append_tag.arg      = args[2]
tags_to_append.arg  = args[3]
append_tags_at.arg  = args[4]
save.at.arg         = args[5]
"
eval(parse(text = arguments.list))

arguments.list = unlist(strsplit(arguments.list, "\n"))
arguments.list = arguments.list[!(arguments.list == "")]

for(n in 1:length(arguments.list)){
  argument = arguments.list[n]
  argument = gsub(pattern=" ", replacement="", x=argument)
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  argument.content = eval(parse(text = argument.name))
  eval(parse(text = argument.content))
  if (!exists(variable.name)){
    stop(sprintf("Argument %s not passed. Stopping ... ", variable.name))
  }
}

for(n in 1:length(seurat.addrs)){
  seurat.addrs[n] = file.path("../../data", seurat.addrs[n])
}
save.at = file.path("../../data", save.at)

if(append_tag){
  if ( length(tags_to_append) != length(seurat.addrs)){
    stop("Number of tags is different from number of data objects. Stopping ... ")
  }
}

library(Seurat)

#######################################################################################################

# init empty list to store the seurat objects
store = list()

# for each .RDS filename read the data and append object to the store lirest
print("loading the datasets")
for (i in 1:length(seurat.addrs)){
  sprintf("Loaded %d out of %d", i, length(seurat.addrs))
  seurat.obj = readRDS(seurat.addrs[i])
  if(append_tag){
    own.tag = tags_to_append[i]
    for(md.index in 1:length(append_tags_at)){
      md.name = append_tags_at[md.index]
      code.line = sprintf("seurat.obj@meta.data$%s = paste(seurat.obj@meta.data$%s, '%s', sep = '_')", md.name, md.name, own.tag)
      eval(parse(text = code.line))
    }
  }
  store[[i]] = seurat.obj
  print(store[[i]])
}

# merge first 2 objects
print("Merging the first 2 datasets")
seurat.obj = MergeSeurat(object1=store[[1]], object2=store[[2]], project="None", min.cells=0, min.genes=0)
# add the rest of the seurat objects
if(length(store) > 2){
  for (j in 3:length(store)){
    sprintf("Adding dataset %d", j)
    seurat.obj = MergeSeurat(object1=seurat.obj, object2=store[[j]], project="None", min.cells=0, min.genes=0)
  }
}
rm(store)

print("normalizing data ... ")
seurat.obj = NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
print("Computing variable genes ... ")
seurat.obj = FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                                dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                x.high.cutoff = 3, y.cutoff = .625)
print("Scaling data ...")
seurat.obj = ScaleData(object=seurat.obj)
print("Saving data ... ")
saveRDS(seurat.obj, save.at)

print("Ended beautifully ... ")
