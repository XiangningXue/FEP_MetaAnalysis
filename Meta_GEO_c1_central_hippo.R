#Meta_GEO_c1_PreDiabetes 
rm(list = ls())
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
hippo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hippocampus")
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
# Prepare for data --------------------------------------------------------
data.contrast = readxl::read_excel(paste0(hippo.dir, "/Full Sample Contrasts_Central Insulin_Hippocampus.xlsx"), skip = 1) 
study.start.rowind = which(!is.na(data.contrast$`Dataset Name`))
study.end.rowind = c(study.start.rowind[-1]-1, nrow(data.contrast))
study.names = data.contrast$`Dataset Name`[study.start.rowind]
study.info = lapply(1:length(study.names), function(s){
  study.info = data.contrast[study.start.rowind[s]:study.end.rowind[s], ]
  info.list = apply(study.info, 2, function(one.col){
    one.col[!is.na(one.col)]
  })
  names(info.list) = colnames(study.info)
  return(info.list)
})
names(study.info) = study.names
# "GSE88723_MEK"      "GSE50873_AMPK_3d"  "GSE50873_AMPK_7d"  "GSE50873_AMPK_14d"
study.info.highlevel = list(GSE88723 = c(1), 
                            GSE50873 = c(2, 3, 4), 
                            GSE63174 = c(5:14),
                            GSE116813 = c(15), 
                            GSE57823 = c(16, 17),
                            GSE64607 = c(18, 19))
library(GEOquery)
library(Biobase)
library(limma)
library(umap)
fdata.gene.col = c("ILMN_Gene", "ILMN_Gene", "ILMN_Gene", "ILMN_Gene", "ILMN_Gene")#only 2 and 4 confirmed
# Do DE for each of the data ----------------------------------------------

# s = 1 -------------------------------------------------------------------
s = 1 #this data set is discard because there is no overlap between the FEP genes
gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset)
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
data.fdata = fData(gset[[idx]])
data.exprs = exprs(gset[[idx]])

# retrieve genes from the gene_assignment column: second item in the // // cell
genes_assignment = data.fdata$ge
genes_assignment2 = lapply(genes_assignment, strsplit, "//")
library(stringr)
genes = sapply(genes_assignment2, function(agene){
  print(paste0(str_trim(agene[[1]][2]), "&&", agene[[1]][2]))
  return(str_trim(agene[[1]][2]))
})
genes = toupper(genes)

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[ex == 0] <- NaN
  ex <- log2(ex) }else{
    ex = ex[!is.na(genes), ]
    genes.filter = genes[!is.na(genes)]
    na.filter = apply(ex, 1, function(a){mean(is.na(a))})
    ex = ex[na.filter<0.5, ]
    genes.filter = genes.filter[na.filter<0.5] 
  }

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", annotation(gset[[idx]]), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", annotation(gset[[idx]]), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

#normalization
ex.norm = ex #no need for normalization
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ex.norm = normalizeBetweenArrays(as.matrix(ex))
# par(mar=c(4,4,2,1))
# title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL6883", " value distribution", sep ="")
# plotDensities(ex.norm, main=title, legend=F)


# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#filter genes that only in the FEP result
sum(rownames(ex2)%in%names(FEP.result)) #2844
ex2.sub = ex2[,c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)),
                         rep(1, length(study.info[[s]]$`treatment samples`))))
colnames(design) = c("Intercept", "trt")
fit = eBayes(lmFit(ex2.sub, design = design))
fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")

write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", study.names[s], ".csv"))


# s = 2, 3, 4 -------------------------------------------------------------------
s = 2 #2, 3, 4 are the same data
  gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")
  
  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])
  
  all(rownames(data.exprs)==rownames(data.fdata))
  which(!toupper(data.fdata$ILMN_Gene)==toupper(data.fdata$Symbol)) #found some Sep genes, but it does not matter because there are no confusing SEP genes in FEP result
  genes = data.fdata$Symbol
  genes[grep("Sep", genes)]
  genes[grep("([0-9]+)-Sep",  genes)]
  genes[grep("([0-9]+)-Sep",  genes)] = paste0("Sep", gsub("([0-9]+)-Sep", "\\1", genes[grep("([0-9]+)-Sep",  genes)]))
  genes = toupper(genes)
  
  ex = data.exprs
  #log2 transformation
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex[ex <= 0] <- NaN
    ex <- log2(ex) 
    }else{
      ex = ex[!is.na(genes), ]
      genes.filter = genes[!is.na(genes)]
      na.filter = apply(ex, 1, function(a){mean(is.na(a))})
      ex = ex[na.filter<0.5|genes%in%names(FEP.result), ]
      genes.filter = genes.filter[na.filter<0.5|genes%in%names(FEP.result)] 
    }
  
  # box-and-whisker plot
  par(mar=c(7,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
  plotDensities(ex, main=title, legend=F)
  
  #normalization -> no need becasue the distribution is identical \
  ex.norm = ex
  # ex <- ex[!duplicated(ex), ]  # remove duplicates
  # ex.norm = normalizeBetweenArrays(as.matrix(ex))
  # par(mar=c(4,4,2,1))
  # title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
  # plotDensities(ex.norm, main=title, legend=F)
  # 
  
  # Filter the duplicated genes ---------------------------------------------
  data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
  ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
    if(nrow(one.gene)==1){
      return(one.gene)
    }else{
      row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
      one.gene = one.gene[which.max(row.cv), ]
      return(one.gene)
    }
  }))
  rownames(ex2) = names(data.exprs.split)
  
  #filter genes that only in the FEP result
  sum(rownames(ex2)%in%names(FEP.result)) #2412
  
#I can just use all three and deal with the multiple testing using limma contrast
#check if the four groups are all different 
length(study.info[[2]]$`control samples`)+length(study.info[[2]]$`treatment samples`)+
  length(study.info[[3]]$`treatment samples`)+length(study.info[[4]]$`treatment samples`) ==
  length(unique(c(study.info[[2]]$`control samples`, study.info[[2]]$`treatment samples`, 
                  study.info[[3]]$`treatment samples`, study.info[[4]]$`treatment samples`)))#TRUE

ex2.sub = ex2[, c(study.info[[2]]$`control samples`, study.info[[2]]$`treatment samples`, 
                          study.info[[3]]$`treatment samples`, study.info[[4]]$`treatment samples`)]
group = factor(c(rep("Control", length(study.info[[2]]$`control samples`)), 
                 rep("AMPK_3d", length(study.info[[2]]$`treatment samples`)), 
                 rep("AMPK_7d", length(study.info[[2]]$`treatment samples`)), 
                 rep("AMPK_14d", length(study.info[[2]]$`treatment samples`))), levels = c("Control", "AMPK_3d", "AMPK_7d", "AMPK_14d"))
design = model.matrix(~group)
fit = eBayes(lmFit(ex2.sub, design = design))
fit.top = topTable(fit, c(2, 3, 4), nrow(fit))

write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", "GSE50873_AMPK_all", ".csv"))

#if doing the correction way 
sep.list = replicate(3, list)
i = 1
for(s in 2:4){
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  sep.list[[i]] = list(res = fit.top, 
                       study = study.names[s])
  i = i+1
  write.csv(fit.top, paste0(hippo.dir, "/all_separate/DE_GEO_", study.names[s], ".csv"))
}

set.seed(15213)
pvalue <- runif(1000, min=0, max=1)
p.list = replicate(4, list)
for(i in 1:length(sep.list)){
  one.res = sep.list[[i]]$res
  ps = one.res$P.Value
  ps[is.na(ps)]=0.5
  p.list[[i]] = gg_qqplot(ps) +
    ggtitle(paste0(sep.list[[i]]$study, ": ", 
                   length(intersect(one.res[, 1], names(FEP.result))), "/", nrow(one.res)))+
    theme_bw(base_size = 12) +
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.15,
      vjust = 1 + 0.15 * 3,
      label = sprintf("λ = %.2f", inflation(ps)),
      size = 4
    ) +
    theme(
      axis.ticks = element_line(size = 0.5),
      panel.grid = element_blank()
      # panel.grid = element_line(size = 0.5, color = "grey80")
    )
}

#correct for the p-value 
genes.sort = sort(rownames(sep.list[[1]]$res))
p.combined = data.frame(p_3d = sep.list[[1]]$res$P.Value[match(genes.sort, rownames(sep.list[[1]]$res))], 
                        p_7d = sep.list[[2]]$res$P.Value[match(genes.sort, rownames(sep.list[[2]]$res))],
                        p_14d = sep.list[[3]]$res$P.Value[match(genes.sort, rownames(sep.list[[3]]$res))])
rownames(p.combined) = genes.sort
p.combined$p.corrected = pbeta(pmin(p.combined$p_3d, p.combined$p_7d, p.combined$p_14d), 1, 3, lower.tail = TRUE)


ps = p.combined$p.corrected
ps[is.na(ps)]=0.5
p.list[[4]] = gg_qqplot(ps) +
  ggtitle(paste0("all"))+
  theme_bw(base_size = 12) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("λ = %.2f", inflation(ps)),
    size = 4
  ) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )

pdf(paste0(hippo.dir, "/all_separate/DE_GEO_", study.names[2], "_qqplot.pdf"), width = 15)
grid.arrange(grobs = p.list, ncol = 2)
dev.off()

all.res = read.csv(paste0(hippo.dir, "/DE_GEO_GSE50873_AMPK_all.csv"))
p.combined$all = all.res$P.Value[match(genes.sort, all.res$X)]
pdf(paste0(hippo.dir, "/all_separate/DE_GEO_", study.names[2], "_pval_compare.pdf"))
plot(-log10(p.combined$all), -log10(p.combined$p.corrected))
curve((x), add = TRUE, col = "green")
dev.off()
plot(p.combined$all, p.combined$p.corrected)
curve((x), add = TRUE, col = "green")

write.csv(p.combined, paste0(hippo.dir, "/all_separate/DE_GEO_GSE50873_AMPK_all.csv"))
# s = 5, 6, ...14 ---------------------------------------------------------
s = 5 
gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset) 
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
data.fdata = fData(gset[[idx]])
data.exprs = exprs(gset[[idx]])

all(rownames(data.exprs)==rownames(data.fdata))
genes = toupper(data.fdata$GENE_SYMBOL)

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[ex <= 0] <- NaN
  ex <- log2(ex) }else{
    ex = ex[!is.na(genes), ]
    genes.filter = genes[!is.na(genes)]
    na.filter = apply(ex, 1, function(a){mean(is.na(a))})
    ex = ex[na.filter<0.5, ]
    genes.filter = genes.filter[na.filter<0.5] 
  }

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

#normalization -> no need becasue the distribution is identical \
#ex.norm = ex
ex.norm = normalizeBetweenArrays(as.matrix(ex))
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(ex.norm, main=title, legend=F)

# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#filter genes that only in the FEP result
sum(rownames(ex2)%in%names(FEP.result)) #3248
for(s in 5:14){
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  
  write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", study.names[s], ".csv"))
}


# s = 15 ------------------------------------------------------------------
s = 15 #2, 3, 4 are the same data
gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset) 
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
data.fdata = fData(gset[[idx]])
data.exprs = exprs(gset[[idx]])

data.exprs = read.table(paste0(hippo.dir, "/GSE116813_ExpectedCounts_genes.out.txt"), header = TRUE, row.names = 1)

library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters = listFilters(ensembl)
filters[grep("NR", filters$description), ]
attributes = listAttributes(ensembl)
attributes[grepl("SYMBOL", toupper(attributes$name)), ]
genes0 = rownames(data.exprs)
genes0 = lapply(genes0, function(a){
  unlist(strsplit(a, ","))
})
genes1 = unlist(genes0)
genesBM1 = getBM(attributes=c('refseq_mrna', "refseq_ncrna", 'mgi_symbol'), 
                filters = c('refseq_mrna'), 
                values = genes1, 
                mart = ensembl)
genesBM2 = getBM(attributes=c('refseq_mrna', "refseq_ncrna", 'mgi_symbol'), 
                 filters = c('refseq_ncrna'), 
                 values = genes1, 
                 mart = ensembl) 
#now take only the first
genes2 = rownames(data.exprs)
genes2 = sapply(genes2, function(a){
  unlist(strsplit(a, ","))[1]
})
match.tab = data.frame(data.gene = toupper(genes2), 
                       symbol1 = toupper(genesBM1[match(genes2, genesBM1$refseq_mrna), "mgi_symbol"]), 
                       symbol2 = toupper(genesBM2[match(genes2, genesBM2$refseq_ncrna), "mgi_symbol"]))
genes = ifelse(!is.na(match.tab$symbol1), match.tab$symbol1, ifelse(!is.na(match.tab$symbol2), match.tab$symbol2, match.tab$data.gene))
names(genes) = rownames(match.tab)

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex.cpm = edgeR::cpm(ex, log = TRUE, prior = 1)
  countCPM=apply(ex.cpm, 1, function(x) mean(x>1))
  ex= ex[countCPM>0.5|(genes%in%names(FEP.result)), ]
  genes.filter = genes[countCPM>0.5|(genes%in%names(FEP.result))]
}else{
  ex = ex[!is.na(genes), ]
  genes.filter = genes[!is.na(genes)]
  na.filter = apply(ex, 1, function(a){mean(is.na(a))})
  ex = ex[na.filter<0.5, ]
  genes.filter = genes.filter[na.filter<0.5]
}

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, sep ="")
boxplot(log2(ex+1), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(log2(ex+1), main=title, legend=F)

#normalization -> no need becasue the distribution is identical \
d = DGEList(counts = ex)
y = voom(d)
ex.norm = y$E
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(ex.norm, main=title, legend=F)

# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#I can not do the analysis becasuse the column names of the actual data does not match what is given. So. I asked Jiwon and now I am going to clean the "GSE116813_CB1R_3wks.xlsx" file
GSE116813 = readxl::read_xlsx(paste0(hippo.dir, "/GSE116813_CB1R_3wks.xlsx"), skip = 1)
colnames(GSE116813) = c("Symbol", "Description", "log2FC", "P.Value", "FDR")
GSE116813$Symbol = toupper(GSE116813$Symbol)
write.csv(GSE116813, paste0(hippo.dir, "/DE_GEO_GSE116813_CB1R_3wks_fromPaper.csv"), row.names = FALSE)

# s = 16, 17 ---------------------------------------------------------
s = 16 #16, 17 are the same data, different samples
gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset) 
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
data.fdata = fData(gset[[idx]])
data.exprs = exprs(gset[[idx]])

all(rownames(data.exprs)==rownames(data.fdata))
genes = toupper(data.fdata$ILMN_Gene)

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[ex == 0] <- NaN
  ex <- log2(ex) }else{
    ex = ex[!is.na(genes), ]
    genes.filter = genes[!is.na(genes)]
    na.filter = apply(ex, 1, function(a){mean(is.na(a))})
    ex = ex[na.filter<0.5, ]
    genes.filter = genes.filter[na.filter<0.5] 
  }

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

#normalization -> no need becasue the distribution is identical \
ex.norm = ex
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ex.norm = normalizeBetweenArrays(as.matrix(ex))
# par(mar=c(4,4,2,1))
# title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
# plotDensities(ex.norm, main=title, legend=F)

# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#filter genes that only in the FEP result
sum(rownames(ex2)%in%names(FEP.result)) #3044
for(s in 16:17){
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  
  write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", study.names[s], ".csv"))
}

# s = 18, 19 ---------------------------------------------------------
s = 18 #16, 17 are the same data, different samples
gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset) 
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
data.fdata = fData(gset[[idx]])
data.exprs = exprs(gset[[idx]])

all(rownames(data.exprs)==rownames(data.fdata))
genes = toupper(data.fdata$ILMN_Gene)

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[ex == 0] <- NaN
  ex <- log2(ex) }else{
    ex = ex[!is.na(genes), ]
    genes.filter = genes[!is.na(genes)]
    na.filter = apply(ex, 1, function(a){mean(is.na(a))})
    ex = ex[na.filter<0.5, ]
    genes.filter = genes.filter[na.filter<0.5] 
  }

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

#normalization -> no need becasue the distribution is identical \
ex.norm = ex
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ex.norm = normalizeBetweenArrays(as.matrix(ex))
# par(mar=c(4,4,2,1))
# title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
# plotDensities(ex.norm, main=title, legend=F)

# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#filter genes that only in the FEP result
sum(rownames(ex2)%in%names(FEP.result)) #2412

for(s in 18:19){
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")

  write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", study.names[s], ".csv"))
}

# New data 20210409--------------------------------------------------------
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
hypo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hypothalamus")
hippo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hippocampus")
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
# Prepare for data --------------------------------------------------------
data.contrast = readxl::read_excel(paste0(dir, "/FEP Metabolic Dysfunction Project", "/Sample contrasts_search update_included studies_Apr0321.xlsx"),  sheet = 2) 
colnames(data.contrast) = c("Dataset Name", "GSE ID", "control samples", "control sample title", "treatment samples", "Treatment sample titles", "Note")
study.start.rowind = which(!is.na(data.contrast$`Dataset Name`))
study.end.rowind = c(study.start.rowind[-1]-1, nrow(data.contrast))
study.names = data.contrast$`Dataset Name`[study.start.rowind]
study.info = lapply(1:length(study.names), function(s){
  study.info = data.contrast[study.start.rowind[s]:study.end.rowind[s], ]
  info.list = apply(study.info, 2, function(one.col){
    one.col[!is.na(one.col)]
  })
  names(info.list) = colnames(study.info)
  return(info.list)
})
names(study.info) = study.names
same.geo = data.contrast$`GSE ID`[!is.na(data.contrast$`GSE ID`)]
study.info.highlevel = list(GSE167264 = c(1), 
                            GSE154434 = c(2))

library(GEOquery)
library(Biobase)
library(limma)
library(biomaRt)

# Do DE for each of the data ----------------------------------------------

# s = 1 -------------------------------------------------------------------
s = 1
gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset) 
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])#there are in total 84 samples
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
# data.fdata = fData(gset[[idx]])
data.exprs = read.table(paste0(hypo.dir, "/GSE167264_Ob_Mus_RNA_seq_counts.txt"), header = TRUE)
colnames(data.exprs) = gsub("(.*)\\.bam", "\\1", colnames(data.exprs))
all(data.pdata$title%in%colnames(data.exprs))#TRUE

#convert emsembl to gene symbol
genes.ensembl = data.exprs$Geneid
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters = listFilters(ensembl)
filters[grep("ENSMUSG", filters$description), ]
attributes = listAttributes(ensembl)
attributes[grepl("SYMBOL", toupper(attributes$name)), ]
genesBM = getBM(attributes=c('ensembl_gene_id', 'mgi_symbol', 'uniprot_gn_symbol'), 
                filters = 'ensembl_gene_id', 
                values = genes.ensembl, 
                mart = ensembl)
data.exprs = data.exprs[data.exprs$Geneid%in%genesBM$ensembl_gene_id, ]
genes = toupper(genesBM[match(data.exprs$Geneid, genesBM$ensembl_gene_id), "mgi_symbol"])
data.exprs = data.exprs[genes!="", ]
genes = genes[genes!=""]

data.exprs = data.exprs[, -c(1:6)]

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex.cpm = edgeR::cpm(ex, log = TRUE, prior = 1)
  countCPM=apply(ex.cpm, 1, function(x) mean(x>1))
  ex= ex[countCPM>0.5|(genes%in%names(FEP.result)), ]
  genes.filter = genes[countCPM>0.5|(genes%in%names(FEP.result))]
}else{
  ex = ex[!is.na(genes), ]
  genes.filter = genes[!is.na(genes)]
  na.filter = apply(ex, 1, function(a){mean(is.na(a))})
  ex = ex[na.filter<0.5, ]
  genes.filter = genes.filter[na.filter<0.5]
}

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, sep ="")
boxplot(log2(ex+1), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(log2(ex+1), main=title, legend=F)

#normalization -> no need becasue the distribution is identical \
#  ex <- ex[!duplicated(ex), ]  # remove duplicates
d = DGEList(counts = ex)
y = voom(d)
ex.norm = y$E
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(ex.norm, main=title, legend=F)


# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#filter genes that only in the FEP result
sum(rownames(ex2)%in%names(FEP.result)) #3154

ex2.sub = ex2[, c(study.info[[1]]$`control sample title`, study.info[[1]]$`Treatment sample titles`)]
group = factor(c(rep("ND", length(study.info[[1]]$`control samples`)), 
                 rep("HFD", length(study.info[[1]]$`treatment samples`))), levels = c("ND", "HFD"))
design = model.matrix(~group)
fit = eBayes(lmFit(ex2.sub, design = design))
fit.top = topTable(fit, c(2), nrow(fit))

write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", study.names[s], ".csv"))

# s = 2 -------------------------------------------------------------------
s = 2
gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
show(gset) 
attr(gset, "names")

#choose GPL
idx <- 1

show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])#there are in total 84 samples
#  table(pData(phenoData(gset[[1]]))[1:24,8])
data.pdata = pData(gset[[idx]])
# data.fdata = fData(gset[[idx]])
data.exprs = readxl::read_xlsx(paste0(hippo.dir, "/GSE154434_FPKM.xlsx"))

sample.ids = colnames(data.exprs)[-1]
genes0 = data.exprs$geneid %>% sapply(function(a){strsplit(a, "\\|")[[1]][2]})
genes = toupper(genes0)

data.exprs = data.exprs[, -1]

ex = data.exprs
#log2 transformation
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex.cpm = edgeR::cpm(ex, log = TRUE, prior = 1)
  countCPM=apply(ex.cpm, 1, function(x) mean(x>1))
  ex= ex[countCPM>0.5|(genes%in%names(FEP.result)), ]
  genes.filter = genes[countCPM>0.5|(genes%in%names(FEP.result))]
}else{
  ex = ex[!is.na(genes), ]
  genes.filter = genes[!is.na(genes)]
  na.filter = apply(ex, 1, function(a){mean(is.na(a))})
  ex = ex[na.filter<0.5, ]
  genes.filter = genes.filter[na.filter<0.5]
}

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, sep ="")
boxplot(log2(ex+1), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(log2(ex+1), main=title, legend=F)

#normalization -> no need becasue the distribution is identical \
#  ex <- ex[!duplicated(ex), ]  # remove duplicates
d = DGEList(counts = ex)
y = voom(d)
ex.norm = y$E
par(mar=c(4,4,2,1))
title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
plotDensities(ex.norm, main=title, legend=F)


# Filter the duplicated genes ---------------------------------------------
data.exprs.split = split(as.data.frame(ex.norm), genes.filter)
ex2 = do.call(rbind.data.frame, lapply(data.exprs.split, function(one.gene){
  if(nrow(one.gene)==1){
    return(one.gene)
  }else{
    row.cv = apply(one.gene, 1, function(a){abs(sd(a, na.rm = TRUE)/mean(a, na.rm = TRUE))})
    one.gene = one.gene[which.max(row.cv), ]
    return(one.gene)
  }
}))
rownames(ex2) = names(data.exprs.split)

#filter genes that only in the FEP result
sum(rownames(ex2)%in%names(FEP.result)) #3084

ex2.sub = ex2[, c(study.info[[s]]$`control sample title`, study.info[[s]]$`Treatment sample titles`)]
group = factor(c(rep("ND", length(study.info[[s]]$`control samples`)), 
                 rep("HFD", length(study.info[[s]]$`treatment samples`))), levels = c("ND", "HFD"))
design = model.matrix(~group)
fit = eBayes(lmFit(ex2.sub, design = design))
fit.top = topTable(fit, c(2), nrow(fit))

write.csv(fit.top, paste0(hippo.dir, "/DE_GEO_", study.names[s], ".csv"))
