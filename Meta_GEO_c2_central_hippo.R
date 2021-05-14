library(dplyr)
library(AWFisher)
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
hippo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hippocampus")
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))

# hippo.data.info = read.csv(paste0(hippo.dir, "/Central_Insulin_Hippocampus_file_dir.csv"))
# hippo.data.info = hippo.data.info[hippo.data.info$dir!="", ]

hippo.files = list.files(hippo.dir)
hippo.DE.file = hippo.files[grepl("DE_GEO", hippo.files)]
hippo.DE.file.name = gsub("DE_GEO_(.*)\\.csv", "\\1", hippo.DE.file)

# # only use the first batch of data
# hippo.DE.file = hippo.DE.file[!(grepl("GSE16724_HFD_Hippo", hippo.DE.file)|grepl("GSE154434_HFD", hippo.DE.file))]
# hippo.DE.file.name = gsub("DE_GEO_(.*)\\.csv", "\\1", hippo.DE.file)

hippo.data.list = lapply(hippo.DE.file, function(dir){
    one.dat = read.csv(paste0(hippo.dir, "/", dir))
    return(one.dat)
})#19 results

names(hippo.data.list) = hippo.DE.file.name
# #replace column name if there is ï
# hippo.data.list = lapply(hippo.data.list, function(one.dat){
#   question.col = grep("ï..", colnames(one.dat))
#   colnames(one.dat)[question.col] =  gsub("ï\\.\\.(.*)", "\\1", colnames(one.dat)[question.col])
#   return(one.dat)
# })
hippo.ID.colnames = c("X")
hippo.gene.list = lapply(1:length(hippo.data.list), function(dat.ind){
  if(sum(grepl("tbl", class(hippo.data.list[[dat.ind]])))){
    gene.list = hippo.data.list[[dat.ind]]%>%pull(hippo.ID.colnames[dat.ind])
  }else{
    gene.list = as.character(hippo.data.list[[dat.ind]][, 1])
  }
})
names(hippo.gene.list) = hippo.DE.file.name

hippo.unique.genes = sort(unique(toupper(unlist(hippo.gene.list))))#length(hippo.unique.genes) #20792
length(hippo.unique.genes) # 49158

#
#Make a table - extract pvals
hippo.ID.colnames = c("X", "Symbol")
hippo.lfc.colnames = c("log2FC", "F", "logFC")
hippo.pval.colnames = c("P.Value")
hippo.lfc.pval.list = lapply(1:length(hippo.data.list), function(dat.ind){
  
  ID.colname = colnames(hippo.data.list[[dat.ind]])[colnames(hippo.data.list[[dat.ind]])%in%hippo.ID.colnames]
  lfc.colname =  colnames(hippo.data.list[[dat.ind]])[colnames(hippo.data.list[[dat.ind]])%in%hippo.lfc.colnames]
  pull.list = hippo.data.list[[dat.ind]][, c(ID.colname, lfc.colname, "P.Value")]
  #}
  colnames(pull.list) = c("Gene.Symbol", lfc.colname,  "P.Value")
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
                          paste0(lfc.colname,"_", names(hippo.data.list)[dat.ind]), 
                          paste0("P.Value_", names(hippo.data.list)[dat.ind]))
  print(dat.ind)
  return(pull.list)
})

hippo.lfc.tab = Reduce(full_join, lapply(hippo.lfc.pval.list, function(one.res){one.res[, -3]}))
hippo.pval.tab = Reduce(full_join, lapply(hippo.lfc.pval.list, function(one.res){one.res[, -2]}))
hippo.lfc.tab = hippo.lfc.tab[order(hippo.lfc.tab$Gene.Symbol), ]
hippo.pval.tab = hippo.pval.tab[order(hippo.pval.tab$Gene.Symbol), ]
sum.value = apply(hippo.lfc.tab[, -1], 1, function(a){sum(!is.na(a))})

# an all-information table 
hippo.lfc.pval.tab = cbind.data.frame(hippo.lfc.tab, hippo.pval.tab[, -1], n.report = sum.value)

#I will filter out only genes in the FEP
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
length(FEP.result)#3727
sum(hippo.pval.tab$Gene.Symbol%in%names(FEP.result)) #3333 
hippo.pval.tab2 = hippo.pval.tab[complete.cases(hippo.pval.tab), ] #5559 -> 7639(in the first batch)
intersect.genes = intersect(hippo.pval.tab2$Gene.Symbol, names(FEP.result)) #1895 -> 2012(in the first batch)
hippo.pval.tab2 = hippo.pval.tab2[hippo.pval.tab2$Gene.Symbol%in%intersect.genes, ]
hippo.pval.tab2 = as.data.frame(hippo.pval.tab2)
rownames(hippo.pval.tab2) = hippo.pval.tab2$Gene.Symbol
hippo.pval.tab2 = hippo.pval.tab2[, -1]

# # Do AW fisher ------------------------------------------------------------
hippo.res = AWFisher_pvalue(as.matrix(hippo.pval.tab2))
hippo.res$gene.names = rownames(hippo.pval.tab2)
# saveRDS(hippo.res, paste0(dir, "/meta_out/meta_p_hippo.rds"))
hippo.res = readRDS(paste0(dir, "/meta_out/meta_p_hippo.rds"))

hippo.res.pvalues = hippo.res$pvalues
names(hippo.res.pvalues) = rownames(hippo.pval.tab2)
FEP.result.sub = FEP.result[intersect.genes]

meta.p.mat = cbind(FEP.result.sub, hippo.res.pvalues)
meta.res =  AWFisher_pvalue(meta.p.mat)
meta.res.q = p.adjust(meta.res$pvalues, "BH")
sum(meta.res$pvalues<0.05) #1323 -> 938 
sum(meta.res.q<0.05) #1173 -> 711

meta.res.tab = data.frame(Genes = intersect.genes,
                          FEP.p = meta.p.mat[, 1],
                          hippo.p = meta.p.mat[, 2],
                          meta.p = meta.res$pvalues,
                          meta.q = meta.res.q,
                          weight11 = apply(meta.res$weights, 1, sum)==2,
                          weight.FEP = meta.res$weights[, 1],
                          weight.T2D = meta.res$weights[, 2])
write.csv(meta.res.tab,
          paste0(dir, "/meta_out/meta_p_FEP&hippo_20210412.csv"))





meta.res.tab = data.frame(Genes = intersect.genes, 
                          FEP.p = FEP.result[intersect.genes], 
                          hippo.p = meta.p.mat[, 2], 
                          meta.p = meta.res$pvalues, 
                          meta.q = meta.res.q, 
                          weight = apply(meta.res$weights, 1, sum)==2)

write.csv(meta.res.tab[meta.res.tab$weight==1, ], 
          paste0(dir, "/meta_out/meta_p_FEP&hippo.csv"))
write.csv(meta.res.tab, 
          paste0(dir, "/meta_out/meta_p_FEP&hippo2.csv"))
#write.csv(meta.res.tab$Genes, 
#          paste0(dir, "/meta_out/meta_p_FEP&hippo_backgroundGenes.csv"))

meta.res.tab = data.frame(Genes = intersect.genes, 
                          FEP.p = FEP.result[intersect.genes], 
                          T2D.p = T2D.res.pvalues[intersect.genes], 
                          meta.p = meta.res$pvalues, 
                          meta.q = meta.res.q, 
                          weight11 = apply(meta.res$weights, 1, sum)==2, 
                          weight.FEP = meta.res$weights[, 1],
                          weight.T2D = meta.res$weights[, 2])

# output meta results only within hippo -----------------------------------
hippo.pval.tab2 = hippo.pval.tab[complete.cases(hippo.pval.tab), ] #3200
rownames(hippo.pval.tab2) = hippo.pval.tab2$Gene.Symbol
hippo.pval.tab2 = hippo.pval.tab2[, -1]
hippo.res = AWFisher_pvalue(as.matrix(hippo.pval.tab2))
hippo.pval.tab2$pval_meta = hippo.res$pvalues
hippo.pval.tab2$qval.meta = p.adjust(hippo.res$pvalues, "BH")
hippo.pval.tab2 = hippo.pval.tab2[order(hippo.pval.tab2$pval_meta), ]
write.csv(hippo.pval.tab2, paste0(dir, "/meta_out/meta_p_hippo.csv"))

#The hippo allinfo output------------
hippo.lfc.pval.tab2 = hippo.lfc.pval.tab
hippo.lfc.pval.tab2 = hippo.lfc.pval.tab2[!(is.na(hippo.lfc.pval.tab2$Gene.Symbol)|hippo.lfc.pval.tab2$Gene.Symbol==""), ]
hippo.lfc.pval.tab2 = hippo.lfc.pval.tab2[apply(hippo.lfc.pval.tab2[, -1], 1, function(a){mean(is.na(a))})!=1, ] #if a gene has no pval or lfc
hippo.lfc.pval.tab2$overlap.gene.FEP = hippo.lfc.pval.tab2$Gene.Symbol%in%intersect.genes
hippo.weight = data.frame(Gene.Symbol = hippo.res$gene.names, hippo.res$weights)
colnames(hippo.weight)[-1] = paste0("weight", "_", names(hippo.data.list))
hippo.lfc.pval.tab2 = left_join(hippo.lfc.pval.tab2, hippo.weight)
hippo.lfc.pval.tab2$meta.pval.hippo = hippo.res$pvalues[match(hippo.lfc.pval.tab2$Gene.Symbol, hippo.res$gene.names)]
hippo.lfc.pval.tab2$meta.pval.FEP = FEP.res[match(hippo.lfc.pval.tab2$Gene.Symbol, names(FEP.res))]
hippo.lfc.pval.tab2$meta.pval.FEP.hippo = meta.res$pvalues[match(hippo.lfc.pval.tab2$Gene.Symbol, rownames(meta.p.mat))]
hippo.lfc.pval.tab2$meta.qval.FEP.hippo = meta.res.q[match(hippo.lfc.pval.tab2$Gene.Symbol, rownames(meta.p.mat))]
meta.weight = data.frame(Gene.Symbol = rownames(meta.p.mat), meta.res$weights)
colnames(meta.weight)[2:3] = c("weight_FEP", "weight_hippo")
meta.weight$weight11 = apply(meta.res$weights, 1, sum)==2
hippo.lfc.pval.tab2 = left_join(hippo.lfc.pval.tab2, meta.weight)

write.csv(hippo.lfc.pval.tab2, paste0(dir, "/meta_out/Meta_p_hippo_allinfo.csv"))

hippo.lfc.pval.tab2 %>%
  filter(overlap.gene.FEP ==1&meta.pval.FEP.hippo<0.05)%>%
  head

#check 
hippo.res$weights[hippo.res$pvalues<0.05, ] %>% colMeans() #the 
hippo.res$weights[hippo.res$pvalues<0.05, ] %>% rowMeans() %>% table
rowmeans0.05 = rowMeans(hippo.res$weights[hippo.res$pvalues<0.05, ])
hippo.res$weights[hippo.res$pvalues<0.05, ][rowmeans0.05<0.06, ]
hippo.res$weights[hippo.res$pvalues<0.05, ][rowmeans0.05>0.06&rowmeans0.05<0.11, ]
hippo.pval.tab2[hippo.res$pvalues<0.05, ][rowmeans0.05<0.06, ]
hippo.pval.tab2[hippo.res$pvalues<0.05, ][rowmeans0.05>0.06&rowmeans0.05<0.11, ]
