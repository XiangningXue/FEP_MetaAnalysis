#Meta_GEO_c1_PreDiabetes 
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
hypo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hypothalamus")
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
# Prepare for data --------------------------------------------------------
data.contrast = readxl::read_excel(paste0(hypo.dir, "/Full Sample Contrasts_Central Insulin_Hypothalamus.xlsx"), skip = 1) 
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
study.info.highlevel = list(GSE104338 = c(1, 2), 
                            GSE73436 = c(3), 
                            GSE104709 = c(4),
                            GSE130597 = c(5), #discarded because it is single cell format. 
                            GSE127056 = c(6, 7),
                            GSE113943 = c(8, 9))

library(GEOquery)
library(Biobase)
library(limma)
library(biomaRt)

fdata.gene.col = c("ILMN_Gene", "ILMN_Gene", "ILMN_Gene", "ILMN_Gene", "ILMN_Gene")#only 2 and 4 confirmed
# Do DE for each of the data ----------------------------------------------
s = 1
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")

  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[1:5,c(1,6,8)])
#  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])
  
# Filter the duplicated genes ---------------------------------------------
  genes_assignment = data.fdata$gene_assignment
  genes_assignment2 = lapply(genes_assignment, strsplit, "//")
  library(stringr)
  genes = sapply(genes_assignment2, function(agene){
    print(paste0(str_trim(agene[[1]][2]), "&&", agene[[1]][2]))
    return(str_trim(agene[[1]][2]))
  })
  genes = toupper(genes)
  
  ex = data.exprs[!is.na(genes), ]
  genes = genes[!is.na(genes)]
  
  #log2 transformation
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    ex[ex == 0] <- NaN
    ex <- log2(ex) }

  # box-and-whisker plot
  par(mar=c(7,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL6883", sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good

  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL6883", " value distribution", sep ="")
  plotDensities(ex, main=title, legend=F)

  #normalization
  ex.filter <- ex[apply(ex, 1, function(x){mean(is.na(x))})<0.5, ]  #remove genes with NA in more than 0.5 samples
  genes.filter <- genes[apply(ex, 1, function(x){mean(is.na(x))})<0.5]
  ex.norm = normalizeBetweenArrays(as.matrix(ex.filter))
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL6883", " value distribution", sep ="")
  plotDensities(ex.norm, main=title, legend=F)

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
  sum(rownames(ex2)%in%names(FEP.result)) #3247

  #combine the two analysis
  length(study.info[[1]]$`control samples`)+length(study.info[[1]]$`treatment samples`)+
    length(study.info[[2]]$`treatment samples`) ==
    length(unique(c(study.info[[1]]$`control samples`, study.info[[1]]$`treatment samples`, 
                    study.info[[2]]$`treatment samples`)))#TRUE
  
  ex2.sub = ex2[, c(study.info[[1]]$`control samples`, study.info[[1]]$`treatment samples`, 
                            study.info[[2]]$`treatment samples`)]
  group = factor(c(rep("Control", length(study.info[[1]]$`control samples`)), 
                   rep("Flaxseed", length(study.info[[1]]$`treatment samples`)), 
                   rep("Safflower", length(study.info[[2]]$`treatment samples`))), levels = c("Control", "Flaxseed", "Safflower"))
  design = model.matrix(~group)
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, c(2, 3), nrow(fit))
  
  write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", "GSE104338_all", ".csv"))

  #try separate test and correct
  sep.list = replicate(2, list)
  i = 1
  for(s in 1:2){
    ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
    design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                             rep(1, length(study.info[[s]]$`treatment samples`))))
    colnames(design) = c("Intercept", "trt")
    fit = eBayes(lmFit(ex2.sub, design = design))
    fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
    sep.list[[i]] = list(res = fit.top, 
                         study = study.names[s])
    i = i+1
    write.csv(fit.top, paste0(hypo.dir, "/all_separate/DE_GEO_", study.names[s], ".csv"))
  }

  set.seed(15213)
  pvalue <- runif(1000, min=0, max=1)
  p.list = replicate(3, list)
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
  p.combined = data.frame(p_Flaxseed = sep.list[[1]]$res$P.Value[match(genes.sort, rownames(sep.list[[1]]$res))], 
                          p_Safflower = sep.list[[2]]$res$P.Value[match(genes.sort, rownames(sep.list[[2]]$res))])
  rownames(p.combined) = genes.sort
  p.combined$p.corrected = pbeta(pmin(p.combined$p_Flaxseed, p.combined$p_Safflower), 1, 2, lower.tail = TRUE)
  
  
  ps = p.combined$p.corrected
  ps[is.na(ps)]=0.5
  p.list[[3]] = gg_qqplot(ps) +
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
  
  pdf(paste0(hypo.dir, "/all_separate/DE_GEO_", study.names[1], "_qqplot.pdf"), width = 15)
  grid.arrange(grobs = p.list, ncol = 3)
  dev.off()
  
  all.res = read.csv(paste0(hypo.dir, "/DE_GEO_GSE104338_all.csv"))
  p.combined$all = all.res$P.Value[match(genes.sort, all.res$X)]
  pdf(paste0(hypo.dir, "/all_separate/DE_GEO_", study.names[1], "_pval_compare.pdf"))
  plot(-log10(p.combined$all), -log10(p.combined$p.corrected))
  curve((x), add = TRUE, col = "green")
  dev.off()
  
  write.csv(p.combined, paste0(hypo.dir, "/all_separate/DE_GEO_GSE104338_all.csv"))
# s = 4 -------------------------------------------------------------------
  s = 3
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")
  # 
  #choose GPL
  idx <- 1
  # 
  # show(pData(phenoData(gset[[idx]]))[1:24,c(1,6,8)])#there are in total 24 samples
  # #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  # data.fdata = fData(gset[[idx]])
  # data.exprs = exprs(gset[[idx]])
  data.exprs = read.table(paste0(hypo.dir, "/GSE73436_data.tsv"), sep = '\t', header = TRUE)
  genes = toupper(data.exprs$gene_short_name)
  data.exprs = data.exprs[, -c(1, 2)]
  
  ex = data.exprs
  #log2 transformation
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex.cpm = edgeR::cpm(ex, log = TRUE, prior = 1)
    countCPM=apply(ex.cpm, 1, function(x) mean(x>1))
    ex= ex[countCPM>0.5|genes%in%names(FEP.result), ]
    genes.filter = genes[countCPM>0.5|genes%in%names(FEP.result)]
    }
  
  # box-and-whisker plot
  par(mar=c(7,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`,  " value distribution", sep ="")
  plotDensities(log2(ex+1), main=title, legend=F)
  
  #normalization -> no need becasue the distribution is identical \
#  ex <- ex[!duplicated(ex), ]  # remove duplicates
  d = DGEList(counts = ex)
  y = voom(d)
  ex.norm = y$E 
  # ex.norm = normalizeBetweenArrays(as.matrix(ex))
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
  sum(rownames(ex2)%in%names(FEP.result)) #3263
  colnames(ex2) = rownames(data.pdata)
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  
  write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", study.names[s], ".csv"))
 

# s=4 ---------------------------------------------------------------------
  s = 4
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset)
  attr(gset, "names")

  #choose GPL
  idx <- 1

  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])
  data.exprs = read.table(paste0(hypo.dir, "/GSE104709_Gene_quant.txt"), sep = '\t', header = TRUE)
  
  # Filter the duplicated genes ---------------------------------------------
  # ensembl <- useMart("ensembl")
  # datasets <- listDatasets(ensembl)
  # datasets[grepl("MUS", toupper(datasets$dataset)), ]
  # ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)

  # ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  # filters = listFilters(ensembl)
  # filters[grep("ENSMUSG", filters$description), ]
  # attributes = listAttributes(ensembl)
  # attributes[grepl("SYMBOL", toupper(attributes$name)), ]
  # genesBM = getBM(attributes=c('ensembl_gene_id', 'mgi_symbol', 'uniprot_gn_symbol'), 
  #       filters = 'ensembl_gene_id', 
  #       values = data.exprs$id, 
  #       mart = ensembl)
  # data.exprs = data.exprs[data.exprs$id%in%genesBM$ensembl_gene_id, ]
  # genes = toupper(genesBM[match(data.exprs$id, genesBM$ensembl_gene_id), "mgi_symbol"])
  genes = toupper(data.exprs$geneSymbol)
  data.exprs = data.exprs[, -1]
  data.exprs = data.exprs[, -22]
  
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
    d = DGEList(counts = ex)
    y = voom(d)
    
  }
  
  # box-and-whisker plot
  par(mar=c(7,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`,  sep ="")
  boxplot(log2(ex+1), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`,  " value distribution", sep ="")
  plotDensities(log2(ex+1), main=title, legend=F)
  
  #normalization -> no need becasue the distribution is identical \
  ex.norm = y$E 
  # ex.norm = normalizeBetweenArrays(as.matrix(ex))
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
  sum(rownames(ex2)%in%names(FEP.result)) #3261
  colnames(ex2) = colnames(ex)
  
  colnames(ex2)[grepl("\\.0([1-9])",colnames(ex2))] =   
    gsub("\\.0([1-9])", "\\.\\1", colnames(ex2)[grepl("\\.0([1-9])",colnames(ex2))])

  colnames(ex2) = rownames(data.pdata)[match(colnames(ex2), make.names(data.pdata$title))]
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  
  write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", study.names[s], ".csv"))
  
# s = 5 -------------------------------------------------------------------
  #discarded because it is single cell format. 

# s = 6 -------------------------------------------------------------------
  s = 6
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset)
  attr(gset, "names")
  
  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])

  # Filter the duplicated genes ---------------------------------------------
  genes_assignment = data.fdata$gene_assignment
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
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
  plotDensities(ex, main=title, legend=F)
  
  #normalization -> no need becasue the distribution is identical \
  #  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ex.norm = ex
  # ex.norm = normalizeBetweenArrays(as.matrix(ex))
  # par(mar=c(4,4,2,1))
  # title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
  # plotDensities(ex.norm, main=title, legend=F)
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
  sum(rownames(ex2)%in%names(FEP.result)) #3198
  
  length(study.info[[6]]$`control samples`)+length(study.info[[6]]$`treatment samples`)+
    length(study.info[[7]]$`treatment samples`) ==
    length(unique(c(study.info[[6]]$`control samples`, study.info[[6]]$`treatment samples`, 
                    study.info[[7]]$`treatment samples`)))#TRUE
  
  ex2.sub = ex2[, c(study.info[[6]]$`control samples`, study.info[[6]]$`treatment samples`, 
                            study.info[[7]]$`treatment samples`)]
  group = factor(c(rep("normal", length(study.info[[6]]$`control samples`)), 
                   rep("HFD4", length(study.info[[6]]$`treatment samples`)), 
                   rep("HFD8", length(study.info[[7]]$`treatment samples`))), levels = c("normal", "HFD4", "HFD8"))
  design = model.matrix(~group)
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, c(2, 3), nrow(fit))
  
  write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", "GSE127056_HFD_all", ".csv"))

  #try separate test and correct
  sep.list = replicate(2, list)
  i = 1
  for(s in 6:7){
    ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
    design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                             rep(1, length(study.info[[s]]$`treatment samples`))))
    colnames(design) = c("Intercept", "trt")
    fit = eBayes(lmFit(ex2.sub, design = design))
    fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
    sep.list[[i]] = list(res = fit.top, 
                         study = study.names[s])
    i = i+1
    write.csv(fit.top, paste0(hypo.dir, "/all_separate/DE_GEO_", study.names[s], ".csv"))
  }
  
  set.seed(15213)
  pvalue <- runif(1000, min=0, max=1)
  p.list = replicate(3, list)
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
  p.combined = data.frame(p_4wk = sep.list[[1]]$res$P.Value[match(genes.sort, rownames(sep.list[[1]]$res))], 
                          p_8wk = sep.list[[2]]$res$P.Value[match(genes.sort, rownames(sep.list[[2]]$res))])
  rownames(p.combined) = genes.sort
  p.combined$p.corrected = pbeta(pmin(p.combined$p_4wk, p.combined$p_8wk), 1, 2, lower.tail = TRUE)
  
  
  ps = p.combined$p.corrected
  ps[is.na(ps)]=0.5
  p.list[[3]] = gg_qqplot(ps) +
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
  
  pdf(paste0(hypo.dir, "/all_separate/DE_GEO_", study.names[6], "_qqplot.pdf"), width = 15)
  grid.arrange(grobs = p.list, ncol = 3)
  dev.off()
  
  all.res = read.csv(paste0(hypo.dir, "/DE_GEO_GSE127056_HFD_all.csv"))
  p.combined$all = all.res$P.Value[match(genes.sort, all.res$X)]
  pdf(paste0(hypo.dir, "/all_separate/DE_GEO_", study.names[6], "_pval_compare.pdf"))
  plot(-log10(p.combined$all), -log10(p.combined$p.corrected))
  curve((x), add = TRUE, col = "green")
  dev.off()
  
  write.csv(p.combined, paste0(hypo.dir, "/all_separate/DE_GEO_GSE127056_HFD_all.csv"))
  
# s = 8 ----------------------------------------------------------------
  s = 8
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset)
  attr(gset, "names")
  
  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])#there are in total 24 samples
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])
  genes = toupper(data.fdata$`External name`)

  ex = data.exprs
  #log2 transformation
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex.cpm = edgeR::cpm(ex, log = TRUE, prior = 1)
    countCPM=apply(ex.cpm, 1, function(x) mean(x>1))
    ex= ex[countCPM>0.5|(rownames(ex)%in%names(FEP.result)), ]
    ex= edgeR::cpm(ex, log = TRUE, prior = 1)
  }else{
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
  #  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ex.norm = ex
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
  sum(rownames(ex2)%in%names(FEP.result)) #2497
  
  for(s in 8:9){
    ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
    group = factor(c(rep("LF", length(study.info[[s]]$`control samples`)), 
                     rep("HF", length(study.info[[s]]$`treatment samples`))), levels = c("LF", "HF"))
    design = model.matrix(~group)
    fit = eBayes(lmFit(ex2.sub, design = design))
    fit.top = topTable(fit, c(2), nrow(fit))
    
    write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", study.names[s], ".csv"))
    
  }
  

# New data 20210409--------------------------------------------------------
  dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
  hypo.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/Central Insulin Datasets/Hypothalamus")
  FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
  # Prepare for data --------------------------------------------------------
  data.contrast = readxl::read_excel(paste0(dir, "/FEP Metabolic Dysfunction Project", "/Sample contrasts_search update_included studies_Apr0321.xlsx"),  sheet = 1) 
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
                              GSE157077 = c(2:7), 
                              GSE145840 = c(8),
                              GSE148641 = c(9:10))
  
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
  data.fdata = fData(gset[[idx]])
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
  
  write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", "GSE167264", ".csv"))
  
  
 # s = 2:7 ----------------------------------------------------------------
  s = 2
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset)
  attr(gset, "names")
  
  #choose GPL
  if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])#there are in total 24 samples
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  # data.fdata = fData(gset[[idx]])
  # data.exprs = exprs(gset[[idx]])
  data.exprs1 = read.table(file = paste0(hypo.dir, "/GSE157077_mouse_scn_control.tsv"), sep = '\t', header = TRUE) 
  data.exprs2 = read.table(file = paste0(hypo.dir, "/GSE157077_mouse_scn_high_fat.tsv"), sep = '\t', header = TRUE) 
  control.time = gsub("ZT_([0-9]+)_REP_.", "\\1", colnames(data.exprs1)[-1])
  hfd.time = gsub("ZT_([0-9]+)_REP_.", "\\1", colnames(data.exprs2)[-1])
  colnames(data.exprs1)[-1] = paste0("ND_", colnames(data.exprs1[-1]))
  colnames(data.exprs2)[-1] = paste0("HFD_", colnames(data.exprs2[-1]))
  
  all(data.exprs1$ID==data.exprs2$ID)
  genes = toupper(data.exprs1$ID)
  length(genes)==length(unique(genes))
  data.exprs = cbind.data.frame(data.exprs1[, -1], data.exprs2[, -1])
  time.point = factor(c(control.time, hfd.time), levels = c(0, 4, 8, 12, 16, 20))
  
  ex = data.exprs
  #log2 transformation
  # log2 transform
  #This is sequencing data, should be count, but now it is not. The lowest number is 0. 
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex.cpm = edgeR::cpm(ex, log = TRUE, prior = 1)
    countCPM=apply(ex.cpm, 1, function(x) mean(x>1))
    ex= ex[countCPM>0.5|genes%in%names(FEP.result), ]
    genes.filter = genes[countCPM>0.5|genes%in%names(FEP.result)]
  }
  
  # box-and-whisker plot
  par(mar=c(7,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", sep ="")
  boxplot(log2(ex+1), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  d = DGEList(counts = ex)
  y = voom(d)
  ex.norm = y$E
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
  plotDensities(ex.norm, main=title, legend=F)
  boxplot(ex.norm, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
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
  sum(rownames(ex2)%in%names(FEP.result)) #3227
  
  group = factor(c(rep("ND", 18), 
                   rep("HFD", 18)), levels = c("ND", "HFD"))
  unique.t = as.character(unique(time.point))
  for(t in 1:length(unique.t)){
    one.time = unique.t[t]
    ex2.sub = ex2[time.point==one.time]
    group.t = group[time.point==one.time]
    design = model.matrix(~group.t)
    
    fit = eBayes(lmFit(ex2.sub, design = design))
    fit.top = topTable(fit, c(2), nrow(fit))
    write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", "GSE157077_HFD_", one.time, "h.csv"))
  }
  
  # s = 8 -------------------------------------------------------------------
  s = 8
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")
  
  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])#there are in total 84 samples
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  # data.fdata = fData(gset[[idx]])
  # data.exprs = exprs(gset[[idx]])

  data.singles = lapply(c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`), function(a){
    a.name = list.files(paste0(hypo.dir, "/GSE145840_RAW"))[grep(a, list.files(paste0(hypo.dir, "/GSE145840_RAW")))]
    one.sample = read.table(paste0(hypo.dir, "/GSE145840_RAW/", a.name))
    colnames(one.sample) = c("id", a)
    return(one.sample)
  })
  gene.overlap = intersect(data.singles[[1]]$id, data.singles[[5]]$id)#already checked that 1-4 have the same genes, and 5-8 have the same
  data.exprs = do.call(cbind.data.frame, lapply(data.singles, function(a){
    b = a[match(gene.overlap, a$id), 2, drop = FALSE]
    rownames(b) = a$id[match(gene.overlap, a$id)]
    return(b)
  }))
  genes = toupper(gene.overlap)
  
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
  # d = DGEList(counts = ex)
  # y = voom(d)
  # ex.norm = y$E
  ex.norm = log2(ex+1)
  ex.norm = normalizeBetweenArrays(ex.norm)
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, " value distribution", sep ="")
  plotDensities(ex.norm, main=title, legend=F, col = c(rep(2, 4), rep(3, 4)))
  boxplot(ex.norm, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  
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
  sum(rownames(ex2)%in%names(FEP.result)) #3201
  
  ex2.sub = ex2[, c(study.info[[8]]$`control samples`, study.info[[8]]$`treatment samples`)]
  group = factor(c(rep("Control", length(study.info[[8]]$`control samples`)), 
                   rep("HFD", length(study.info[[8]]$`treatment samples`))), levels = c("Control", "HFD"))
  design = model.matrix(~group)
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, c(2), nrow(fit))

  # sort(fit.top[names(FEP.result),  "P.Value"])
 
  # ps = fit.top[names(FEP.result),  "P.Value"][complete.cases(fit.top[names(FEP.result),  "P.Value"])]
  # gg_qqplot(ps) +
  #   ggtitle(paste0(result.names[i], ": ",
  #                  length(intersect(one.res[, 1], names(FEP.result))), "/", nrow(one.res)))+
  #   theme_bw(base_size = 12) +
  #   annotate(
  #     geom = "text",
  #     x = -Inf,
  #     y = Inf,
  #     hjust = -0.15,
  #     vjust = 1 + 0.15 * 3,
  #     label = sprintf("λ = %.2f", inflation(ps)),
  #     size = 4
  #   ) +
  #   theme(
  #     axis.ticks = element_line(size = 0.5),
  #     panel.grid = element_blank()
  #     # panel.grid = element_line(size = 0.5, color = "grey80")
  #   )
  
  write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", "GSE145840_HFD", ".csv"))
  # fit.top = read.csv(paste0(hypo.dir, "/DE_GEO_GSE104338_all.csv"))
  # qqPlot(fit.top$P.Value, "unif")
  
  # s = 9:10 -------------------------------------------------------------------
  s = 9
  gset <- getGEO(study.info[[s]]$`GSE ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")
  
  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[,c(1,6,8)])#there are in total 84 samples
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  # data.fdata = fData(gset[[idx]])
  # data.exprs = exprs(gset[[idx]])
  
  data.exprs1 = read.table(paste0(hypo.dir, "/GSE148641_NC_ARC_VMH_counts.txt"), header = TRUE)
  data.exprs2 = read.table(paste0(hypo.dir, "/GSE148641_hypothalamus_reverb_HFD_circadian_RNA_seq.txt"), header = TRUE)
  gene.overlap = intersect(data.exprs1$Geneid, data.exprs2$Geneid)
  data.exprs1 = data.exprs1[match(gene.overlap, data.exprs1$Geneid), ]
  data.exprs2 = data.exprs2[match(gene.overlap, data.exprs2$Geneid), ]

  data.exprs = cbind.data.frame(data.exprs1[, -c(1:6)], data.exprs2[, -c(1:6)])
  data.exprs = data.exprs[, c(study.info[[9]]$`control sample title`, 
                              study.info[[9]]$`Treatment sample titles`,
                              study.info[[10]]$`control sample title`, 
                              study.info[[10]]$`Treatment sample titles`)]
  genes = toupper(gene.overlap)
  
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
  sum(rownames(ex2)%in%names(FEP.result)) #3227
  
  study.name = character(10)
  study.name[9:10] = c("ARC", "VMH")
  for(s in 9:10){
    ex2.sub = ex2[, c(study.info[[s]]$`control sample title`, study.info[[s]]$`Treatment sample titles`)]
    group = factor(c(rep("Control", length(study.info[[s]]$`control sample title`)), 
                     rep("HFD", length(study.info[[s]]$`Treatment sample titles`))), levels = c("Control", "HFD"))
    design = model.matrix(~group)
    fit = eBayes(lmFit(ex2.sub, design = design))
    fit.top = topTable(fit, c(2), nrow(fit))
    
    write.csv(fit.top, paste0(hypo.dir, "/DE_GEO_", "GSE148641_", study.name[s], ".csv"))
  }
  
  
  