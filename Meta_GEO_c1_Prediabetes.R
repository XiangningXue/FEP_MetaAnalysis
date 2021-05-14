#Meta_GEO_c1_PreDiabetes 
dir = "/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Logan/FEP_diabetes" 
T2D.dir = paste0(dir, "/FEP Metabolic Dysfunction Project/T2D Datasets")
FEP.result = readRDS(paste0(dir, "/meta_out/Meta_p_FEP.rds"))
# Prepare for data --------------------------------------------------------
data.contrast = readxl::read_excel(paste0(T2D.dir, "/Full Sample Contrasts_T2D_mRNA.xlsx")) 
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

library(GEOquery)
library(Biobase)
library(limma)
library(umap)
fdata.gene.col = c("ILMN_Gene", "ILMN_Gene", "ILMN_Gene", "ILMN_Gene", "ILMN_Gene")#only 2 and 4 confirmed
# Do DE for each of the data ----------------------------------------------
s = 2
  gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")

  #choose GPL
  if (length(gset) > 1) idx <- grep("GPL6883", attr(gset, "names")) else idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[1:24,c(1,6,8)])#there are in total 24 samples
#  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])
  
# Filter the duplicated genes ---------------------------------------------
  genes = data.fdata[rownames(data.exprs), fdata.gene.col[s]] 

  ex = data.exprs
  #log2 transformation
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
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
  ex.filter <- ex[apply(ex, 1, function(x){mean(is.na(x))})<0.5|(genes%in%names(FEP.result)), ]  #remove genes with NA in more than 0.5 samples
  genes.filter <- genes[apply(ex, 1, function(x){mean(is.na(x))})<0.5|(genes%in%names(FEP.result))]
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
  sum(rownames(ex2)%in%names(FEP.result)) #3204
  ex2.sub = ex2[,c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  
  write.csv(fit.top, paste0(T2D.dir, "/DE_GEO_", study.names[s], ".csv"))
  

# s = 4 -------------------------------------------------------------------
  s = 4
  gset <- getGEO(study.info[[s]]$`GEO Accession ID`, GSEMatrix =TRUE, getGPL=TRUE)
  show(gset) 
  attr(gset, "names")
  
  #choose GPL
  idx <- 1
  
  show(pData(phenoData(gset[[idx]]))[1:24,c(1,6,8)])
  #  table(pData(phenoData(gset[[1]]))[1:24,8])
  data.pdata = pData(gset[[idx]])
  data.fdata = fData(gset[[idx]])
  data.exprs = exprs(gset[[idx]])
  
  # Filter the duplicated genes ---------------------------------------------
  genes = data.fdata[rownames(data.exprs), fdata.gene.col[s]] 
  
  ex = data.exprs
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
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)#distribution looks good
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
  plotDensities(ex, main=title, legend=F)
  
  #normalization -> no need becasue the distribution is identical \
  ex.filter <- ex[apply(ex, 1, function(x){mean(is.na(x))})<0.5|(genes%in%names(FEP.result)), ]  #remove genes with NA in more than 0.5 samples
  genes.filter <- genes[apply(ex, 1, function(x){mean(is.na(x))})<0.5|(genes%in%names(FEP.result))]
  ex.norm = ex.filter
  # ex <- ex[!duplicated(ex), ]  # remove duplicates
  # ex.norm = normalizeBetweenArrays(as.matrix(ex))
  # par(mar=c(4,4,2,1))
  # title <- paste (study.info[[s]]$`GEO Accession ID`, "/", "GPL10558", " value distribution", sep ="")
  # plotDensities(ex.norm, main=title, legend=F)
  # 
  
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
  sum(rownames(ex2)%in%names(FEP.result)) #3727
  ex2.sub = ex2[, c(study.info[[s]]$`control samples`, study.info[[s]]$`treatment samples`)]
  design = model.matrix(~c(rep(0, length(study.info[[s]]$`control samples`)), 
                           rep(1, length(study.info[[s]]$`treatment samples`))))
  colnames(design) = c("Intercept", "trt")
  fit = eBayes(lmFit(ex2.sub, design = design))
  fit.top = topTable(fit, 2, nrow(fit), sort.by = "P")
  
  write.csv(fit.top, paste0(T2D.dir, "/DE_GEO_", study.names[s], ".csv"))
  
  