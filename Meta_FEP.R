#macbook
# libraries and directories -----------------------------------------------
library(dplyr)
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
FEP.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/FEP Transcriptomic Datasets")


# Do DE analysis for the raw data -----------------------------------------
#the study is Mantere_WB_Individual Expression.xlsx
dat.Mantere = readxl::read_excel(paste0(FEP.dir, "/Mantere_WB_Individual Expression.xlsx"), skip = 2)
dat.clinical = dat.Mantere[, 1:3]
n.control = nrow(dat.clinical)-4
dat.clinical$Disease = factor(c(rep(0, n.control), rep(1, 4)))
dat.clinical$Sex = factor(dat.clinical$Sex)
dat.exp = t(dat.Mantere[, -c(1:3)])
genes.Mantere = rownames(dat.exp)
dat.exp = apply(dat.exp, 2, as.numeric)
dat.exp = log2(dat.exp)
rownames(dat.exp) = genes.Mantere
colnames(dat.exp) = dat.clinical$Sample

library(limma)
design = model.matrix(~Disease+Sex, data = dat.clinical)
colnames(design)[1] = "Intercept"
fit = eBayes(lmFit(dat.exp, design))
top.fit = topTable(fit, coef = 2, n = nrow(dat.exp), sort.by = "P")
write.csv(top.fit, paste0(FEP.dir, "/Mantere_WB_DE.csv"))

#Kumara
dat.Kumara = readxl::read_excel(paste0(FEP.dir, "/Kumarasinghe-GEX-raw-clean.xlsx"))
genes.Kumara = dat.Kumara$SYMBOL
cols.select.Kumara = colnames(dat.Kumara)[grepl("AVG_Signal", colnames(dat.Kumara))]
cols.select.Kumara2 = cols.select.Kumara[!grepl("after", cols.select.Kumara) ]
dat.Kumara2 = dat.Kumara[, cols.select.Kumara2]
#check dup.genes
table(table(genes.Kumara))
group = c(rep("control", 11), rep("trt", 10))
boxplot(dat.Kumara2, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
ex.norm = normalizeBetweenArrays(as.matrix(dat.Kumara2))
boxplot(ex.norm, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
library(limma)
design = model.matrix(~group)
fit = eBayes(lmFit(ex.norm, design = design))
fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
fit.top[c("135", "136"), ]
fit = eBayes(lmFit(dat.Kumara2, design = design))
fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
fit.top[c("135", "136"), ]
mean(as.numeric(dat.Kumara2[135, 12:21]))/mean(as.numeric(dat.Kumara2[135, 1:11]))
mean(as.numeric(dat.Kumara2[136, 12:21]))/mean(as.numeric(dat.Kumara2[136, 1:11]))

mean(as.numeric(dat.Kumara2[135,  1:11]))/mean(as.numeric(dat.Kumara2[135, 12:21]))
mean(as.numeric(dat.Kumara2[136,  1:11]))/mean(as.numeric(dat.Kumara2[136, 12:21]))

mean(as.numeric(dat.Kumara2[135,  1:11]))/mean(as.numeric(dat.Kumara2[135, 12:21]))
log2((2^mean(log2(as.numeric(dat.Kumara2[135,  1:11]))))/(2^mean(log2(as.numeric(dat.Kumara2[135, 12:21])), na.rm = TRUE)))
log2((2^mean(log2(as.numeric(dat.Kumara2[136,  1:11]))))/(2^mean(log2(as.numeric(dat.Kumara2[136, 12:21])))))
(2^mean(log2(as.numeric(dat.Kumara2[136,  1:11]))))/(2^mean(log2(as.numeric(dat.Kumara2[136, 12:21]))))

mean(as.numeric(dat.Kumara2[229,  1:11]))/mean(as.numeric(dat.Kumara2[229, 12:21]))
mean(as.numeric(dat.Kumara2[230,  1:11]))/mean(as.numeric(dat.Kumara2[230, 12:21]))
(2^mean(log2(as.numeric(dat.Kumara2[230,  1:11]))))/(2^mean(log2(as.numeric(dat.Kumara2[230, 12:21]))))

#sam r
# library(samr)
# samr.data = list(x=as.matrix(dat.Kumara2),y=factor(group), geneid= 1:nrow(dat.Kumara2), genenames= dat.Kumara$SYMBOL, logged2=F)
# samr.obj = samr(samr.data, resp.type="Two class unpaired", nperms=1000, random.seed=15217)
# delta.table = samr.compute.delta.table(samr.obj)
# head(delta.table[delta.table[, "median FDR"]<0.05, ])
# delta.choose = delta.table[delta.table[, "median FDR"]<0.05, "delta"][1]
# siggenes.table = samr.compute.siggenes.table(samr.obj, del=delta.choose, samr.data, delta.table)
# 
# samr.data2 = list(x=as.matrix(dat.Kumara2[1:1000, ]),y=factor(group), geneid= 1:1000, genenames= dat.Kumara$SYMBOL[1:1000], logged2=T)
# samr.obj2 = samr(samr.data2, resp.type="Two class unpaired", nperms=1000, random.seed=15217)
# delta.table = samr.compute.delta.table(samr.obj2)
# head(delta.table[delta.table[, "median FDR"]<0.05, ])
# delta.choose = delta.table[delta.table[, "median FDR"]<0.05, "delta"][1]
# siggenes.table2 = samr.compute.siggenes.table(samr.obj2, del=delta.choose, samr.data, delta.table)


# 
# #decide to use DE from the paper. 
# DE.original.Kumara = readxl::read_xlsx(paste0(FEP.dir, "/Kumarasinghe-GEX-DE.xlsx"), skip = 2)
# sum(is.na(DE.original.Kumara$`Illumnia ID`))#2
# DE.original.Kumara[is.na(DE.original.Kumara$`Illumnia ID`), ]#it is the separation line
# DE.original.Kumara = DE.original.Kumara[!is.na(DE.original.Kumara$`Illumnia ID`), ]
# DE.original.Kumara$`Fold Change (SZ/CTL)` = as.numeric(DE.original.Kumara$`Fold Change (SZ/CTL)`)
# lfc.na.which = which(is.na(DE.original.Kumara$`Fold Change (SZ/CTL)`))
# DE.original.Kumara[lfc.na.which, ]#should not include
# DE.original.Kumara = DE.original.Kumara[!is.na(DE.original.Kumara$`Fold Change (SZ/CTL)`), ]
# DE.original.Kumara$`p-value*` = as.numeric(DE.original.Kumara$`p-value*`)
# write.csv(DE.original.Kumara, paste0(FEP.dir, "/Kumarasinghe-GEX-DE2.csv"))
# dup.genes.Kumara = names(table(DE.original.Kumara$`Gene Symbol`))[table(DE.original.Kumara$`Gene Symbol`)>1]
# DE.original.Kumara[DE.original.Kumara$`Gene Symbol`%in%dup.genes.Kumara, ]

# Duplicated genes in study 4 ---------------------------------------------
study4.DE = FEP.data.list[[4]]
study4.DE = study4.DE[!is.na(study4.DE$SYMBOL), ]
genes.dup5 = names(table(study4.DE$SYMBOL))[table(study4.DE$SYMBOL)==5]
study4.DE[study4.DE$SYMBOL==genes.dup5, ]

# Make a gene * sutdy matrix that contains pvals --------------------------
#read prepared data
FEP.data.info = read.csv(paste0(FEP.dir, "/FEP_file_dir.csv"))
FEP.data.info = FEP.data.info[FEP.data.info$dir!="", ]
FEP.data.list = lapply(FEP.data.info$dir, function(dir){
  if(grepl("csv", dir)){
    one.dat = read.csv(paste0(FEP.dir, "/", dir))
  }else{
    one.dat = readxl::read_excel(paste0(FEP.dir, "/", dir))
  }
})
names(FEP.data.list) = FEP.data.info$Author.year
# #replace column name if there is ï
# FEP.data.list = lapply(FEP.data.list, function(one.dat){
#   question.col = grep("ï..", colnames(one.dat))
#   colnames(one.dat)[question.col] =  gsub("ï\\.\\.(.*)", "\\1", colnames(one.dat)[question.col])
#   return(one.dat)
# })
FEP.ID.colnames = c("TargetID", "Gene.Symbol", "X", "SYMBOL", "Gene", 
                    "Gene", "Gene", "Gene name", "Gene.Symbol")
FEP.gene.list = lapply(1:length(FEP.data.list), function(dat.ind){
    if(sum(grepl("tbl", class(FEP.data.list[[dat.ind]])))){
      gene.list = FEP.data.list[[dat.ind]]%>%pull(FEP.ID.colnames[dat.ind])
    }else{
      gene.list = as.character(FEP.data.list[[dat.ind]][, FEP.ID.colnames[dat.ind]])
    }
})
names(FEP.gene.list) = FEP.data.info$Author.year

FEP.unique.genes = FEP.gene.list %>% 
  unlist %>%
  toupper %>%
  unique %>%
  sort
FEP.unique.genes = FEP.unique.genes[-1]#length(FEP.unique.genes) 23244
#the first one is deleted because it is ""

#after using the DE file for Kumarasinghe2013 the number is 23243
#check what genes have lower case names
lapply(FEP.gene.list, function(one.list){
  one.list[grepl("[a-z]", one.list)]
})
#By the way, look for Sep genes and march
lapply(FEP.gene.list, function(one.list){
  one.list[grepl("Sep", one.list)|grepl("Mar", one.list)|grepl("SEP", one.list)|grepl("MAR", one.list)]
}) #it seem that there is no confusion, only Sainz 2013 has some of these genes

# Check for duplicated genes ----------------------------------------------
FEP.ID.colnames = c("TargetID", "Gene.Symbol", "X", "SYMBOL", "Gene", 
                    "Gene", "Gene", "Gene name", "Gene.Symbol")
lapply(FEP.gene.list, function(a){table(table(a))}) #there are some duplicated genes, need to choose one

#Make a table - extract pvals and lfc
FEP.ID.colnames = c("TargetID", "Gene.Symbol", "X", "SYMBOL", "Gene", 
                    "Gene", "Gene", "Gene name", "Gene.Symbol")
FEP.data.list[[7]]$logFC = log2(FEP.data.list[[7]]$FC)
FEP.lfc.colnames = c("logFC", "log2FoldChange", "logFC", "logFC", "statistic", "statistic", "logFC", "Log2FC", "Fold.Change..SZ.CTL.")
#for study 9. It is fold change, with direction of change added as sign. 
#by the way, what is this statistic?
FEP.pval.colnames = c("P.Value", "pval", "P.Value", "P.Value", "P.value", "P.value", "P value", "P-value", "p.value.")
FEP.lfc.pval.list = lapply(1:length(FEP.pval.colnames), function(dat.ind){
  # if(sum(grepl("tbl", class(FEP.data.list[[dat.ind]])))){
  #   pull.list = FEP.data.list[[dat.ind]]%>%
  #     pull(FEP.ID.colnames[dat.ind], FEP.lfc.colnames[dat.ind], FEP.pval.colnames[dat.ind])
  # }else{
    pull.list = FEP.data.list[[dat.ind]][, c(FEP.ID.colnames[dat.ind], FEP.lfc.colnames[dat.ind], FEP.pval.colnames[dat.ind])]
  #}
  colnames(pull.list) = c("Gene.Symbol", "logFC", "P.Value")
  pull.list$Gene.Symbol = toupper(pull.list$Gene.Symbol)
  pull.list$logFC = as.numeric(pull.list$logFC)
  pull.list$P.Value = as.numeric(pull.list$P.Value)
  pull.list = pull.list[!is.na(pull.list$P.Value), ]
  colnames(pull.list) = c("Gene.Symbol", paste0("logFC_", dat.ind), paste0("P.Value_", dat.ind))
  return(pull.list)
})

FEP.lfc.tab = Reduce(full_join, lapply(FEP.lfc.pval.list, function(one.res){one.res[, -3]}))
FEP.pval.tab = Reduce(full_join, lapply(FEP.lfc.pval.list, function(one.res){one.res[, -2]}))
FEP.lfc.tab = FEP.lfc.tab[order(FEP.lfc.tab$Gene.Symbol), ]
FEP.pval.tab = FEP.pval.tab[order(FEP.pval.tab$Gene.Symbol), ]
sum.value = apply(FEP.lfc.tab[, -1], 1, function(a){sum(!is.na(a))})
FEP.lfc.tab[sum.value==7, ] #a little hard to tell... 

#write out this table with all information
FEP.lfc.pval.tab = cbind.data.frame(FEP.lfc.tab, FEP.pval.tab[, -1], n.report = sum.value)
change.tab = data.frame(ind = 1:9, study.name = names(FEP.data.list), statistics.name = FEP.lfc.colnames)
change.tab$stat_study = paste0(change.tab$statistics.name, "_", change.tab$study.name)
change.tab$pval_study = paste0("P.Value", "_", change.tab$study.name)
colnames(FEP.lfc.pval.tab)[2:19] = c(change.tab$stat_stud, change.tab$pval_study)
use.ind = c(1, 2, 4, 5, 6, 9)

# Do AW-fisher ------------------------------------------------------------
library(AWFisher)
#res = AWFisher_pvalue(FEP.pval.tab)
#R corrupts with the above command
#try two columns first
#res = AWFisher_pvalue(FEP.pval.tab[, 1:2])
#can AWFisher handle NA? 
FEP.pval.tab[sum.value==7, ]
FEP.pval.tab7 = FEP.pval.tab[sum.value==7, c(2, 3, 5, 10)]
FEP.pval.tab7 = FEP.pval.tab[sum.value==7, c(2, 3)]
res = AWFisher_pvalue(FEP.pval.tab7)


qvalue <- p.adjust(res$pvalue, "BH") 
meta_p_q = data.frame(pvalue = res$pvalues, qvalue = qvalue, row.names = row.names(pmatrix))
sigIndex = which(res$weights[,1]==1&res$weights[,2]==1)
meta_p_q = meta_p_q[sigIndex,]
meta_p_q_sort = meta_p_q[order(meta_p_q$qvalue,decreasing = F),]

# trunP: find genes that overlaps in the "all"  ---------------------------
all.ind = c(1, 2, 4)
overlap.gene = Reduce(intersect, lapply(FEP.lfc.pval.list[all.ind], function(g){g$Gene.Symbol}))
length(overlap.gene)#3727

#make pmatrix
FEP.pval.tab2 = FEP.pval.tab[FEP.pval.tab$Gene.Symbol%in%overlap.gene, ]
use.ind = c(1, 2, 4, 5, 6, 9)
FEP.pval.tab2 = FEP.pval.tab2[, c(1, use.ind+1)]
#deal with duplicated gene issue. 
sum(duplicated(FEP.pval.tab2$Gene.Symbol))#737 genes
UniqueGenes.FEP.pval.tab2 = unique(FEP.pval.tab2$Gene.Symbol)
FEP.pval.tab2.listByGene = split(FEP.pval.tab2, FEP.pval.tab2$Gene.Symbol)
FEP.pval.tab2.listKeepOne = lapply(FEP.pval.tab2.listByGene, function(a){
  n = dim(a)[1]
  if(n == 1){
    return(a)
  }else if(n>1){
    mean.by.row = apply(a, 1, function(aa){mean(as.numeric(aa[-1]), na.rm = TRUE)})
    keep = which.min(mean.by.row)
    a2 = a[keep, ]
    return(a2)
  }
})
FEP.pval.tab3 = do.call(rbind.data.frame,FEP.pval.tab2.listKeepOne)
FEP.pval.tab3 = as.data.frame(FEP.pval.tab3)
rownames(FEP.pval.tab3) = FEP.pval.tab3$Gene.Symbol
FEP.pval.tab3 = FEP.pval.tab3[, -1]
FEP.pval.threshold = c(0, 0, 0, 0.01, 0.01, 0.05)


# Mean Imputation -------------------------------------------------
### Generate permutation matrix ########
# i = 1
# one.pvec = FEP.pval.tab3[i, ]
# a = one.pvec[!is.na(one.pvec)]
# alpha = FEP.pval.threshold[is.na(one.pvec)]
mean.Fisher <- function(a,alpha){
  
  ### a: vector of observed p values
###  ### b: vector of unobserved p values
  ### alpha: threshold of unobserved p values
  
  ### length of observed p values ####
  n.1 <- length(a)
  
  ### length of censored p values ####
  n.2 <- length(alpha)
  
  ### define the test statistic t #####
  
  b <- (1+alpha)/2
  x <- c(a,b)
  x <- ifelse(x==0,1e-20,x)
  t <- sum(-2*log(x))
  
  ### define a vector of integers from 0 to n.2 ####
  
  nn <- 0:(n.2)
  nn.grid = expand.grid(replicate(n.2, list(c(0, 1))))
  
  ### define the binomial density at alpha #####
  
  coef.binom <- apply(nn.grid, 1, function(a){prod((alpha^a)*((1-alpha)^(1-a)))})
  
  ### define several constants #####

  c.1 <- 2*log((1+alpha)/2)
  c.2 <- 2*log((1+alpha)/alpha)
  tt.test <- apply(nn.grid, 1, function(a){t+sum(c.1)-sum(a*c.2)})
  
  ### compute chi-squre terms ###
  
  test.stat <- pchisq(tt.test, df=2*n.1,lower.tail=FALSE)
  
  ### compute the cdf #######
  
  t.cdf <- sum((test.stat)*(coef.binom))
  return(t.cdf)
}

#####1. Fisher's method ######################################################

p.mean.Fisher <- apply(FEP.pval.tab3,1,function(one.pvec){
  a = one.pvec[!is.na(one.pvec)]
  alpha = FEP.pval.threshold[is.na(one.pvec)]
  mean.Fisher(a,alpha)
})
q.mean.Fisher <- p.adjust(p.mean.Fisher,method="BH")
sum(q.mean.Fisher<0.05) ## 149
sum(q.mean.Fisher<0.01) ## 40
sum(q.mean.Fisher<0.001) ## 9
sum(q.mean.Fisher<0.0001) ## 3
sum(q.mean.Fisher<1e-5) ## 2
sum(q.mean.Fisher<1e-6) ## 2
sum(q.mean.Fisher<1e-7) ## 1

saveRDS(p.mean.Fisher, paste0(dir, "/meta_out/Meta_p_FEP.rds"))

FEP.lfc.pval.tab2 = FEP.lfc.pval.tab[, c(1, use.ind+1, use.ind+10, 20)]
FEP.lfc.pval.tab2 = FEP.lfc.pval.tab2[!is.na(FEP.lfc.pval.tab2$Gene.Symbol), ]
FEP.lfc.pval.tab2 = FEP.lfc.pval.tab2[apply(FEP.lfc.pval.tab2[, -1], 1, function(a){mean(is.na(a))})!=1, ]
FEP.lfc.pval.tab2$overlap.gene = FEP.lfc.pval.tab2$Gene.Symbol%in%overlap.gene
FEP.lfc.pval.tab2.dup = FEP.lfc.pval.tab2[FEP.lfc.pval.tab2$overlap.gene, ]

UniqueGenes.FEP.pval.tab2 = unique(FEP.lfc.pval.tab2.dup$Gene.Symbol)
FEP.pval.tab2.listByGene2 = split(FEP.lfc.pval.tab2.dup, FEP.lfc.pval.tab2.dup$Gene.Symbol)
FEP.pval.tab2.listKeepOne2 = lapply(FEP.pval.tab2.listByGene2, function(a){
  n = dim(a)[1]
  if(n == 1){
    return(a)
  }else if(n>1){
    mean.by.row = apply(a, 1, function(aa){mean(as.numeric(aa[7:12]), na.rm = TRUE)})
    keep = which.min(mean.by.row)
    a2 = a
    a2[1, ] = a[keep, ]
    a2[-1, ] = -99
    return(a2)
  }
})
FEP.lfc.pval.tab2.dup = do.call(rbind.data.frame,FEP.pval.tab2.listKeepOne2)
FEP.lfc.pval.tab2[FEP.lfc.pval.tab2$overlap.gene, ] =FEP.lfc.pval.tab2.dup
#FEP.lfc.pval.tab2[FEP.lfc.pval.tab2$overlap.gene!=0, ]
FEP.lfc.pval.tab2 = FEP.lfc.pval.tab2[FEP.lfc.pval.tab2$Gene.Symbol!=-99, ]
FEP.lfc.pval.tab2$meta.pval = p.mean.Fisher[match(FEP.lfc.pval.tab2$Gene.Symbol, names(p.mean.Fisher))]

write.csv(FEP.lfc.pval.tab2, paste0(dir, "/meta_out/Meta_p_FEP_allinfo.csv"))

FEP.lfc.pval.tab2 %>%
  filter(overlap.gene ==1&meta.pval<0.05)%>%
  head
