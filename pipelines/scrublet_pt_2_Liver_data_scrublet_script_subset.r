
Seurat=readRDS('./liver.RDS')

Seurat@meta.data$fullmeta=names(Seurat@ident)

table(Seurat@meta.data$fullmeta)

library(stringr)
output=str_split(Seurat@meta.data$fullmeta, "_")

short_code=NULL
for( i in 1:length(output)){
  short_code[i]=paste(output[[i]][4],"_",output[[i]][5],sep="")
}

output=str_split(Seurat@meta.data$fullmeta, "_([^_]*)$")
output[1]
for( i in 1:length(output)){
  output[[i]] = output[[i]][-2]
}


Seurat@meta.data$fullmeta = unlist(output)
table(Seurat@meta.data$fullmeta)

table(Seurat@meta.data$fetal.ids)
table(Seurat@meta.data$sort.ids)
table(Seurat@meta.data$tissue)
table(Seurat@meta.data$stages)
table(Seurat@meta.data$sample.type)
table(Seurat@meta.data$gender)
table(Seurat@meta.data$AnnatomicalPart)
table(Seurat@meta.data$lanes)


VlnPlot(Seurat,'RPS4Y1')
VlnPlot(Seurat,'XIST')

#### Load in doublet output
scrub1 = read.csv('./scrublet-scores/all.csv',header=T,row.names=1)

table(rownames(scrub1) %in% short_code)

scrub1=scrub1[rownames(scrub1) %in% short_code,]

match_up = match(rownames(scrub1),short_code)

scrub1=scrub1[order(match_up),]

Seurat@meta.data$scrublet_score = scrub1$scrublet_score
Seurat@meta.data$scrublet_cluster_score = scrub1$scrublet_cluster_score
Seurat@meta.data$bh_pval = scrub1$bh_pval

table(Seurat@meta.data$bh_pval>0.05)

VlnPlot(Seurat, features.plot = "nGene")
min(Seurat@meta.data$nGene)
Seurat@meta.data$percent.mito


saveRDS(Seurat,"./liver_post_scrub.RDS")

###Optional###
Convert(Seurat, 'anndata', X.slot = "data", raw.slot = "raw.data", filename = './Desktop/Haniffa_Group_Work/Paper - Yolk Sac/Post_liver_paper_work/Liver_for_YS.h5ad')





