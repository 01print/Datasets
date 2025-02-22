############
# 2021-09-04/05 
# small room in 2f // from 21-09-09, in 2#,104-2,XJTU-hc, / 9-10/13 //
# the datasets might needed in the future
# glioma
# Pancheng Wu
###########

## GSE16011
GSE16011 {
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE16011", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ## 
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
} 

## GSE23806  GPL570 
GSE23806 {
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE23806", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ## 
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
}

## GSE43378  GPL570
GSE43378 {
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE43378", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ## 
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
}

## GSE12907  GPL96
GSE43378 {
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE12907", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ## 
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
}

## GSE4412 // GPL96   GPL97
GSE4412 {
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE4412", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset1 <- gset$`GSE4412-GPL96_series_matrix.txt.gz` 
  gset2 <- gset$`GSE4412-GPL97_series_matrix.txt.gz`
  
  ## 表达矩阵 GPL96处理
  expr_GPL96{
    exprSet1 = exprs(gset1)
    a <- exprSet1
    exprSet1 <- normalizeBetweenArrays(exprSet1)  ## 标准化一番
    exprSet1 <- as.data.frame(exprSet1) ## 转为数据框
    exprSet1 <- log2(exprSet1)
    gene_id1 <- rownames(exprSet1)
    exprSet1 <- cbind("probe_id" = gene_id1, exprSet1)
    exprSet1$probe_id <- as.character(exprSet1$probe_id)
    
    ## 
    load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL96.Rda")
    exprSet1 <- exprSet1%>%
      inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
      select(-probe_id) %>% #去掉多余信息
      select(symbol, everything()) %>% #重新排列，
      mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
      arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
      distinct(symbol,.keep_all = T) %>% # symbol留下第一个
      select(-rowMean) %>% #反向选择去除rowMean这一列
      tibble::column_to_rownames(colnames(.)[1])
    
    save(exprSet1, file = "exprSet_GSE4412_GPL96.Rda")  ## 85个
    
  }
  ## 表达矩阵 GPL97处理
  expr_GPL97{
    exprSet2 = exprs(gset2)
    b <- exprSet2
    exprSet2 <- normalizeBetweenArrays(exprSet2)  ## 标准化一番
    exprSet2 <- as.data.frame(exprSet2) ## 转为数据框
    exprSet2 <- log2(exprSet2)
    gene_id2 <- rownames(exprSet2)
    exprSet2 <- cbind("probe_id" = gene_id2, exprSet2)
    exprSet2$probe_id <- as.character(exprSet2$probe_id)
    
    ## 制作probe2symbol, GPL97
    GPL97.17394 <- read.delim("C:/Users/Pancheng/Downloads/GPL97-17394.txt", comment.char="#")
    probe2symbol_GPL97 <- GPL97.17394[,c(1,11)]
    names(probe2symbol_GPL97) <- c("probe_id","symbol")
    save(probe2symbol_GPL97, file = "probe2symbol_GPL97.Rda")
    
    load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol_GPL97.Rda")
    exprSet2 <- exprSet2%>%
      inner_join(probe2symbol_GPL97,by="probe_id") %>% #合并探针的信息
      select(-probe_id) %>% #去掉多余信息
      select(symbol, everything()) %>% #重新排列，
      mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
      arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
      distinct(symbol,.keep_all = T) %>% # symbol留下第一个
      select(-rowMean) %>% #反向选择去除rowMean这一列
      tibble::column_to_rownames(colnames(.)[1])
    
    save(exprSet2, file = "exprSet_GSE4412_GPL97.Rda")  ## 85个, ## 11350 obs
    
  }
  
  ## 临床信息处理
  pdata = pData(gset1)
  pdata1 <- pdata[,c(2,43,45,46,49,50,52)]
  expr1 <- as.data.frame(t(exprSet1))
  expr1 <- expr1[,-1]
  GSM <- rownames(expr1)
  expr1 <- cbind("geo_accession" = GSM,expr1)
  expr1_clinic <- merge(expr1, pdata1, by = "geo_accession")
  save(expr1_clinic, file = "expr_clinic_GSE4412_GPL96.Rda")
  
  #
  pdata <- pData(gset2)
  pdata2 <- pdata[,c(2,43,45,46,49,50,52)]
  expr2 <- as.data.frame(t(exprSet2))
  expr2 <- expr2[,-1]
  GSM <- rownames(expr2)
  expr2 <- cbind("geo_accession" = GSM,expr2)
  expr2_clinic <- merge(expr2, pdata2, by = "geo_accession")
  save(expr2_clinic, file = "expr_clinic_GSE4412_GPL97.Rda") ##似乎样本一样，测序平台不一样
}
