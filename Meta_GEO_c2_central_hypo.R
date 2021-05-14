library(dplyr)
library(AWFisher)
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
hypo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hypothalamus")
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))

# hypo.data.info = read.csv(paste0(hypo.dir, "/Central_Insulin_hypocampus_file_dir.csv"))
# hypo.data.info = hypo.data.info[hypo.data.info$dir!="", ]

hypo.files = list.files(hypo.dir)
hypo.DE.file = hypo.files[grepl("DE_GEO", hypo.files)]
hypo.DE.file.name = gsub("DE_GEO_(.*)\\.csv", "\\1", hypo.DE.file)

# #if only use the first batch
# data.contrast = readxl::read_excel(paste0(hypo.dir, "/Full Sample Contrasts_Central Insulin_Hypothalamus.xlsx"), skip = 1) 
# study.start.rowind = which(!is.na(data.contrast$`Dataset Name`))
# study.end.rowind = c(study.start.rowind[-1]-1, nrow(data.contrast))
# study.names = data.contrast$`Dataset Name`[study.start.rowind]
# study.info = lapply(1:length(study.names), function(s){
#   study.info = data.contrast[study.start.rowind[s]:study.end.rowind[s], ]
#   info.list = apply(study.info, 2, function(one.col){
#     one.col[!is.na(one.col)]
#   })
#   names(info.list) = colnames(study.info)
#   return(info.list)
# })
# names(study.info) = study.names
# same.geo = unique(data.contrast$`GSE ID`[!is.na(data.contrast$`GSE ID`)])
# hypo.DE.file.gse = gsub("(GSE[0-9]+)_.*", "\\1", hypo.DE.file.name)
# hypo.DE.file = hypo.DE.file[hypo.DE.file.gse%in%same.geo]
# hypo.DE.file.name = gsub("DE_GEO_(.*)\\.csv", "\\1", hypo.DE.file)


hypo.data.list = lapply(hypo.DE.file, function(dir){
    one.dat = read.csv(paste0(hypo.dir, "/", dir), row.names = 1)
    one.dat$X = rownames(one.dat)
    return(one.dat)
})#13 results

names(hypo.data.list) = hypo.DE.file.name
# #replace column name if there is ï
# hypo.data.list = lapply(hypo.data.list, function(one.dat){
#   question.col = grep("ï..", colnames(one.dat))
#   colnames(one.dat)[question.col] =  gsub("ï\\.\\.(.*)", "\\1", colnames(one.dat)[question.col])
#   return(one.dat)
# })
hypo.ID.colnames = c("X")
hypo.gene.list = lapply(1:length(hypo.data.list), function(dat.ind){
  if(sum(grepl("tbl", class(hypo.data.list[[dat.ind]])))){
    gene.list = hypo.data.list[[dat.ind]]%>%pull(hypo.ID.colnames[dat.ind])
  }else{
    gene.list = as.character(hypo.data.list[[dat.ind]]$X)
  }
})
names(hypo.gene.list) = hypo.DE.file.name

hypo.unique.genes = sort(unique(toupper(unlist(hypo.gene.list))))#length(hypo.unique.genes) #20792
length(hypo.unique.genes) # 18118->77737 (2021-04-18)

#
#Make a table - extract pvals
hypo.ID.colnames = c("X")
hypo.lfc.colnames = c("F", "logFC")
hypo.pval.colnames = c("P.Value")
hypo.lfc.pval.list = lapply(1:length(hypo.data.list), function(dat.ind){
  # if(sum(grepl("tbl", class(hypo.data.list[[dat.ind]])))){
  #   pull.list = hypo.data.list[[dat.ind]]%>%
  #     pull(hypo.ID.colnames[dat.ind], hypo.lfc.colnames[dat.ind], hypo.pval.colnames[dat.ind])
  # }else{
  
  lfc.colname =  colnames(hypo.data.list[[dat.ind]])[colnames(hypo.data.list[[dat.ind]])%in%hypo.lfc.colnames]
  
  pull.list = hypo.data.list[[dat.ind]][, c("X",lfc.colname, "P.Value")]
  #}
  colnames(pull.list) = c("Gene.Symbol", lfc.colname, "P.Value")
  pull.list$Gene.Symbol = toupper(pull.list$Gene.Symbol)
  #pull.list$F = as.numeric(pull.list$F)
  pull.list$P.Value = as.numeric(pull.list$P.Value)
  pull.list = pull.list[!is.na(pull.list$P.Value), ]
  
  #keep only one record: the most significant ones
  # pull.list.split = split(pull.list, pull.list$Gene.Symbol)
  # pull.list = do.call(rbind.data.frame, lapply(pull.list.split, function(one.gene.list){
  #   keep.one = one.gene.list[which.min(one.gene.list$P.Value), ]
  # }))
  
  colnames(pull.list) = c("Gene.Symbol",
                          paste0(lfc.colname,"_", names(hypo.data.list)[dat.ind]), 
                          paste0("P.Value_", names(hypo.data.list)[dat.ind]))
  print(dat.ind)
  return(pull.list)
})

hypo.lfc.tab = Reduce(full_join, lapply(hypo.lfc.pval.list, function(one.res){one.res[, -3]}))
hypo.pval.tab = Reduce(full_join, lapply(hypo.lfc.pval.list, function(one.res){one.res[, -2]}))
hypo.lfc.tab = hypo.lfc.tab[order(hypo.lfc.tab$Gene.Symbol), ]
hypo.pval.tab = hypo.pval.tab[order(hypo.pval.tab$Gene.Symbol), ]
sum.value = apply(hypo.lfc.tab[, -1], 1, function(a){sum(!is.na(a))})

# an all-information table 
hypo.lfc.pval.tab = cbind.data.frame(hypo.lfc.tab, hypo.pval.tab[, -1], n.report = sum.value)

#I will filter out only genes in the FEP
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
length(FEP.result)#3727
sum(hypo.pval.tab$Gene.Symbol%in%names(FEP.result)) #3361
hypo.pval.tab2 = hypo.pval.tab[complete.cases(hypo.pval.tab), ] #3200 -> 2794 (add second batch)
intersect.genes = intersect(hypo.pval.tab2$Gene.Symbol, names(FEP.result)) #2338
hypo.pval.tab2 = hypo.pval.tab2[hypo.pval.tab2$Gene.Symbol%in%intersect.genes, ]
hypo.pval.tab2 = as.data.frame(hypo.pval.tab2)
rownames(hypo.pval.tab2) = hypo.pval.tab2$Gene.Symbol
hypo.pval.tab2 = hypo.pval.tab2[, -1]

# # Do AW fisher ------------------------------------------------------------
hypo.res = AWFisher_pvalue(as.matrix(hypo.pval.tab2))
hypo.res$gene.names = rownames(hypo.pval.tab2)
saveRDS(hypo.res, paste0(dir, "/meta_out/meta_p_hypo.rds"))

hypo.res = readRDS(paste0(dir, "/meta_out/meta_p_hypo.rds"))

hypo.res.pvalues = hypo.res$pvalues
names(hypo.res.pvalues) = rownames(hypo.pval.tab2)
FEP.result.sub = FEP.result[intersect.genes]

meta.p.mat = cbind(FEP.result.sub, hypo.res.pvalues)
meta.res =  AWFisher_pvalue(meta.p.mat)
meta.res.q = p.adjust(meta.res$pvalues, "BH")
sum(meta.res$pvalues<0.05) #579 -> 458 
sum(meta.res.q<0.05) #116 > 80

meta.res.tab = data.frame(Genes = intersect.genes,
                          FEP.p = meta.p.mat[, 1],
                          hypo.p = meta.p.mat[, 2],
                          meta.p = meta.res$pvalues,
                          meta.q = meta.res.q,
                          weight11 = apply(meta.res$weights, 1, sum)==2,
                          weight.FEP = meta.res$weights[, 1],
                          weight.T2D = meta.res$weights[, 2])
write.csv(meta.res.tab,
          paste0(dir, "/meta_out/meta_p_FEP&hypo_20210412.csv"))

meta.res.tab = data.frame(Genes = intersect.genes, 
                          FEP.p = FEP.result[intersect.genes], 
                          hypo.p = meta.p.mat[, 2], 
                          meta.p = meta.res$pvalues, 
                          meta.q = meta.res.q, 
                          weight = apply(meta.res$weights, 1, sum)==2)

write.csv(meta.res.tab[meta.res.tab$weight==1, ], 
          paste0(dir, "/meta_out/meta_p_FEP&hypo.csv"))
write.csv(meta.res.tab, 
          paste0(dir, "/meta_out/meta_p_FEP&hypo2.csv"))
#write.csv(meta.res.tab$Genes, 
#          paste0(dir, "/meta_out/meta_p_FEP&hypo_backgroundGenes.csv"))

# output meta results only within hypo -----------------------------------
hypo.pval.tab2 = hypo.pval.tab[complete.cases(hypo.pval.tab), ] #3200
rownames(hypo.pval.tab2) = hypo.pval.tab2$Gene.Symbol
hypo.pval.tab2 = hypo.pval.tab2[, -1]
hypo.res = AWFisher_pvalue(as.matrix(hypo.pval.tab2))
hypo.pval.tab2$pval_meta = hypo.res$pvalues
hypo.pval.tab2$qval.meta = p.adjust(hypo.res$pvalues, "BH")
hypo.pval.tab2 = hypo.pval.tab2[order(hypo.pval.tab2$pval_meta), ]
write.csv(hypo.pval.tab2, paste0(dir, "/meta_out/meta_p_hypo.csv"))

#The hypo allinfo output------------
hypo.lfc.pval.tab2 = hypo.lfc.pval.tab
hypo.lfc.pval.tab2 = hypo.lfc.pval.tab2[!(is.na(hypo.lfc.pval.tab2$Gene.Symbol)|hypo.lfc.pval.tab2$Gene.Symbol==""), ]
hypo.lfc.pval.tab2 = hypo.lfc.pval.tab2[apply(hypo.lfc.pval.tab2[, -1], 1, function(a){mean(is.na(a))})!=1, ] #if a gene has no pval or lfc
hypo.lfc.pval.tab2$overlap.gene.FEP = hypo.lfc.pval.tab2$Gene.Symbol%in%intersect.genes
hypo.weight = data.frame(Gene.Symbol = hypo.res$gene.names, hypo.res$weights)
colnames(hypo.weight)[-1] = paste0("weight", "_", names(hypo.data.list))
hypo.lfc.pval.tab2 = left_join(hypo.lfc.pval.tab2, hypo.weight)
hypo.lfc.pval.tab2$meta.pval.hypo = hypo.res$pvalues[match(hypo.lfc.pval.tab2$Gene.Symbol, hypo.res$gene.names)]
hypo.lfc.pval.tab2$meta.pval.FEP = FEP.res[match(hypo.lfc.pval.tab2$Gene.Symbol, names(FEP.res))]
hypo.lfc.pval.tab2$meta.pval.FEP.hypo = meta.res$pvalues[match(hypo.lfc.pval.tab2$Gene.Symbol, rownames(meta.p.mat))]
hypo.lfc.pval.tab2$meta.qval.FEP.hypo = meta.res.q[match(hypo.lfc.pval.tab2$Gene.Symbol, rownames(meta.p.mat))]
meta.weight = data.frame(Gene.Symbol = rownames(meta.p.mat), meta.res$weights)
colnames(meta.weight)[2:3] = c("weight_FEP", "weight_hypo")
meta.weight$weight11 = apply(meta.res$weights, 1, sum)==2
hypo.lfc.pval.tab2 = left_join(hypo.lfc.pval.tab2, meta.weight)

write.csv(hypo.lfc.pval.tab2, paste0(dir, "/meta_out/Meta_p_hypo_allinfo.csv"))

hypo.lfc.pval.tab2 %>%
  filter(overlap.gene.FEP ==1&meta.pval.FEP.hypo<0.05)%>%
  head

#check 
hypo.res$weights[hypo.res$pvalues<0.05, ] %>% colMeans() #the 
hypo.res$weights[hypo.res$pvalues<0.05, ] %>% rowMeans() %>% table
rowmeans0.05 = rowMeans(hypo.res$weights[hypo.res$pvalues<0.05, ])
hypo.res$weights[hypo.res$pvalues<0.05, ][rowmeans0.05<0.08, ]
hypo.res$weights[hypo.res$pvalues<0.05, ][rowmeans0.05>0.08&rowmeans0.05<0.16, ]
hypo.pval.tab2[hypo.res$pvalues<0.05, ][rowmeans0.05<0.06, ]
hypo.pval.tab2[hypo.res$pvalues<0.05, ][rowmeans0.05>0.06&rowmeans0.05<0.11, ]


