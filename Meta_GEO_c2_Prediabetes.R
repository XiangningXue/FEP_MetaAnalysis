library(dplyr)
library(AWFisher)
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
T2D.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/T2D Datasets")

T2D.data.info = read.csv(paste0(T2D.dir, "/T2D_file_dir.csv"))
T2D.data.info = T2D.data.info[T2D.data.info$dir!="", ]
T2D.data.info = T2D.data.info[T2D.data.info$Patient.Type=="Prediabetic", ]

T2D.data.list = lapply(T2D.data.info$dir, function(dir){
  if(grepl("csv", dir)){
    one.dat = read.csv(paste0(T2D.dir, "/", dir))
  }else{
    one.dat = readxl::read_excel(paste0(T2D.dir, "/", dir))
    one.dat = one.dat[complete.cases(one.dat),]
  }
})
names(T2D.data.list) = T2D.data.info$study.name
# #replace column name if there is ï
# T2D.data.list = lapply(T2D.data.list, function(one.dat){
#   question.col = grep("ï..", colnames(one.dat))
#   colnames(one.dat)[question.col] =  gsub("ï\\.\\.(.*)", "\\1", colnames(one.dat)[question.col])
#   return(one.dat)
# })
T2D.ID.colnames = c("X","X")
T2D.gene.list = lapply(1:length(T2D.data.list), function(dat.ind){
  if(sum(grepl("tbl", class(T2D.data.list[[dat.ind]])))){
    gene.list = T2D.data.list[[dat.ind]]%>%pull(T2D.ID.colnames[dat.ind])
  }else{
    gene.list = as.character(T2D.data.list[[dat.ind]][, T2D.ID.colnames[dat.ind]])
  }
})
names(T2D.gene.list) = T2D.data.info$study.name

T2D.unique.genes = sort(unique(toupper(unlist(T2D.gene.list))))#length(T2D.unique.genes) #20792
length(T2D.unique.genes) # 36477 because already filtered

#
#Make a table - extract pvals
T2D.ID.colnames = c("X","X")
T2D.lfc.colnames = c("logFC", "logFC")
T2D.pval.colnames = c("P.Value", "P.Value")
T2D.lfc.pval.list = lapply(1:length(T2D.pval.colnames), function(dat.ind){
  # if(sum(grepl("tbl", class(T2D.data.list[[dat.ind]])))){
  #   pull.list = T2D.data.list[[dat.ind]]%>%
  #     pull(T2D.ID.colnames[dat.ind], T2D.lfc.colnames[dat.ind], T2D.pval.colnames[dat.ind])
  # }else{
  pull.list = T2D.data.list[[dat.ind]][, c(T2D.ID.colnames[dat.ind], T2D.lfc.colnames[dat.ind], T2D.pval.colnames[dat.ind])]
  #}
  colnames(pull.list) = c("Gene.Symbol", "logFC", "P.Value")
  pull.list$Gene.Symbol = toupper(pull.list$Gene.Symbol)
  pull.list$logFC = as.numeric(pull.list$logFC)
  pull.list$P.Value = as.numeric(pull.list$P.Value)
  pull.list = pull.list[!is.na(pull.list$P.Value), ]
  
  colnames(pull.list) = c("Gene.Symbol", paste0("logFC_", dat.ind), paste0("P.Value_", dat.ind))
  return(pull.list)
})

T2D.lfc.tab = Reduce(full_join, lapply(T2D.lfc.pval.list, function(one.res){one.res[, -3]}))
T2D.pval.tab = Reduce(full_join, lapply(T2D.lfc.pval.list, function(one.res){one.res[, -2]}))
T2D.lfc.tab = T2D.lfc.tab[order(T2D.lfc.tab$Gene.Symbol), ]
T2D.pval.tab = T2D.pval.tab[order(T2D.pval.tab$Gene.Symbol), ]
sum.value = apply(T2D.lfc.tab[, -1], 1, function(a){sum(!is.na(a))})

# an all-information table 
T2D.lfc.pval.tab = cbind.data.frame(T2D.lfc.tab, T2D.pval.tab[, -1], n.report = sum.value)
change.tab = data.frame(ind = 1:2, study.name = names(T2D.data.list), statistics.name = T2D.lfc.colnames)
change.tab$stat_study = paste0(change.tab$statistics.name, "_", change.tab$study.name)
change.tab$pval_study = paste0("P.Value", "_", change.tab$study.name)
colnames(T2D.lfc.pval.tab)[2:5] = c(change.tab$stat_stud, change.tab$pval_study)

#I will filter out only genes in the FEP
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
length(FEP.result)#3727
sum(T2D.pval.tab$Gene.Symbol%in%names(FEP.result)) 
T2D.pval.tab2 = T2D.pval.tab[complete.cases(T2D.pval.tab), ] #3200
intersect.genes = intersect(T2D.pval.tab2$Gene.Symbol, names(FEP.result)) #3171
T2D.pval.tab2 = T2D.pval.tab2[T2D.pval.tab2$Gene.Symbol%in%intersect.genes, ]
T2D.pval.tab2 = as.data.frame(T2D.pval.tab2)
rownames(T2D.pval.tab2) = T2D.pval.tab2$Gene.Symbol
T2D.pval.tab2 = T2D.pval.tab2[, -1]
 
# Do AW fisher ------------------------------------------------------------
T2D.res = AWFisher_pvalue(as.matrix(T2D.pval.tab2))
T2D.res$gene.names = rownames(T2D.pval.tab2)
saveRDS(T2D.res, paste0(dir, "/meta_out/meta_p_T2D.rds"))

#
T2D.res = readRDS(paste0(dir, "/meta_out/meta_p_T2D.rds"))
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))

#
T2D.res.pvalues = T2D.res$pvalues
names(T2D.res.pvalues) = T2D.res$gene.names
FEP.result.sub = FEP.result[names(T2D.res.pvalues)]

meta.p.mat = cbind(FEP.result.sub, T2D.res.pvalues)
meta.res =  AWFisher_pvalue(meta.p.mat)
meta.res.q = p.adjust(meta.res$pvalues, "BH")
sum(meta.res$pvalues<0.05) #782 -> 826
sum(meta.res.q<0.05) #175 -> 183
#result updated after I changed the order of the pipeline: the variance is calculated on the normalized and transformed data

# meta.res.tab = data.frame(Genes = intersect.genes, 
#                           FEP.p = FEP.result[intersect.genes], 
#                           T2D.p = T2D.res.pvalues[intersect.genes], 
#                           meta.p = meta.res$pvalues, 
#                           meta.q = meta.res.q, 
#                           weight = apply(meta.res$weights, 1, sum)==2)
# 
# write.csv(meta.res.tab[meta.res.tab$weight==1, ], 
#           paste0(dir, "/meta_out/meta_p_FEP&T2D.csv"))
# write.csv(meta.res.tab, 
#           paste0(dir, "/meta_out/meta_p_FEP&T2D2.csv"))
# write.csv(meta.res.tab$Genes, 
#           paste0(dir, "/meta_out/meta_p_FEP&T2D_backgroundGenes.csv"))
# 
# meta.res.tab = data.frame(Genes = intersect.genes, 
#                           FEP.p = FEP.result[intersect.genes], 
#                           T2D.p = T2D.res.pvalues[intersect.genes], 
#                           meta.p = meta.res$pvalues, 
#                           meta.q = meta.res.q, 
#                           weight11 = apply(meta.res$weights, 1, sum)==2, 
#                           weight.FEP = meta.res$weights[, 1],
#                           weight.T2D = meta.res$weights[, 2])


meta.res.tab = data.frame(Genes = intersect.genes,
                          FEP.p = meta.p.mat[, 1],
                          T2D.p = meta.p.mat[, 2],
                          meta.p = meta.res$pvalues,
                          meta.q = meta.res.q,
                          weight11 = apply(meta.res$weights, 1, sum)==2,
                          weight.FEP = meta.res$weights[, 1],
                          weight.T2D = meta.res$weights[, 2])
write.csv(meta.res.tab,
          paste0(dir, "/meta_out/meta_p_FEP&T2D_20210412.csv"))

#add module colors also.
hdir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes"
WGCNA.dir1.4 = paste0(hdir, "/WGCNA_FEPmeta/WGCNA_module_all/power4")
module.dir = WGCNA.dir1.4
module = read.csv(paste0(module.dir, "/Module_all_samples_cut5.csv"), row.names = 1)
meta.res.tab$module = module$module[match(meta.res.tab$Genes, module$genes)]

#meta.res.tab$module[is.na(meta.res.tab$module)]
write.csv(meta.res.tab,
          paste0(dir, "/meta_out/meta_p_FEP&T2D_20210412_module.csv"))


T2D.lfc.pval.tab2 = T2D.lfc.pval.tab
T2D.lfc.pval.tab2 = T2D.lfc.pval.tab2[!is.na(T2D.lfc.pval.tab2$Gene.Symbol), ]
T2D.lfc.pval.tab2 = T2D.lfc.pval.tab2[apply(T2D.lfc.pval.tab2[, -1], 1, function(a){mean(is.na(a))})!=1, ] #if a gene has no pval or lfc
T2D.lfc.pval.tab2$overlap.gene.FEP = T2D.lfc.pval.tab2$Gene.Symbol%in%intersect.genes
T2D.weight = data.frame(Gene.Symbol = T2D.res$gene.names, T2D.res$weights)
colnames(T2D.weight)[2:3] = paste0("weight", "_", names(T2D.data.list))
T2D.lfc.pval.tab2 = left_join(T2D.lfc.pval.tab2, T2D.weight)
T2D.lfc.pval.tab2$meta.pval.T2D = T2D.res$pvalues[match(T2D.lfc.pval.tab2$Gene.Symbol, T2D.res$gene.names)]
T2D.lfc.pval.tab2$meta.pval.FEP = FEP.res[match(T2D.lfc.pval.tab2$Gene.Symbol, names(FEP.res))]
T2D.lfc.pval.tab2$meta.pval.FEP.T2D = meta.res$pvalues[match(T2D.lfc.pval.tab2$Gene.Symbol, rownames(meta.p.mat))]
T2D.lfc.pval.tab2$meta.qval.FEP.T2D = meta.res.q[match(T2D.lfc.pval.tab2$Gene.Symbol, rownames(meta.p.mat))]
meta.weight = data.frame(Gene.Symbol = rownames(meta.p.mat), meta.res$weights)
colnames(meta.weight)[2:3] = c("weight_FEP", "weight_T2D")
meta.weight$weight11 = apply(meta.res$weights, 1, sum)==2
T2D.lfc.pval.tab2 = left_join(T2D.lfc.pval.tab2, meta.weight)
  
write.csv(T2D.lfc.pval.tab2, paste0(dir, "/meta_out/Meta_p_T2D_allinfo.csv"))

T2D.lfc.pval.tab2 %>%
  filter(overlap.gene ==1&meta.pval<0.05)%>%
  head
