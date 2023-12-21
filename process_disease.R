
rm(list=ls())
library(dplyr)
library(data.table)
library(reshape2)
library(reticulate)
use_python('/home/huzixin/anaconda3/envs/test/bin/python')
# file_path <- '/mnt/Shared_Data/ewas/disease_methylation_v1.txt'
# load('/mnt/Shared_Data/ewas/sample_disease.RData')

chunk_size <- 10000 # 可根据需要调整大小
chunk_start <- 1

results <- list()
datas <- list()

repeat {
  print(paste(chunk_start,Sys.time()))
  data <- fread(input = file_path, skip = chunk_start - 1, nrows = chunk_size, header = TRUE)
  
  if (nrow(data) == 0) {
    break
  }
  
  if(chunk_start==1){
    map <- t(data[1:3,-1])
    map <- data.table(
      sample = rownames(map),
      disease = map[,1],
      tissue = map[,2],
      type = map[,3]
    ) %>%
      mutate(key=paste(type,disease,ifelse(tissue=='whole blood','blood','tissue'),sep='_'))
    map$age <- sample_disease$age

    data <- data[-1:-3,]
  }
  
  X <- data[,-1] %>% as.matrix
  rownames(X) <- data[[1]]
  X <- X[rowMeans(is.na(X))<=0.2,,drop=F]
  Xdnames <- dimnames(X)
  X <- apply(X,2,as.numeric)
  dimnames(X) <- Xdnames
  # X <- apply(X,2,function(x){ifelse(is.na(x),median(x),x)})
  X[is.na(X)] <- 0
  X <- log(X+1e-6)
  
  agep <- apply(X,1,function(x){cor.test(x,map$age,use='pairwise')$p.value})
  ager <- cor(t(X),map$age,use='pairwise')
  ager <- data.table(site=rownames(X),p=agep,r=as.vector(ager)) %>% mutate(padj=p.adjust(p,method='fdr'))
  
  Xmean <- apply(X,1,function(x){
    tapply(x,map$key,mean)
  }) %>%
    melt() %>%
    select(key=1,site=2,mean=value)
  Xsd <- apply(X,1,function(x){
    tapply(x,map$key,sd)
  }) %>%
    melt() %>%
    select(key=1,site=2,sd=value)
  Xn <- apply(X,1,function(x){
    tapply(x,map$key,length)
  }) %>%
    melt() %>%
    select(key=1,site=2,n=value)
  out <- merge(Xmean,Xsd) %>% merge(Xn)
  out <- do.call(rbind,strsplit(paste(out$key),'_')) %>%
    as.data.table() %>%
    select(case=1,disease=2,tissue=3) %>%
    data.table(out[,-1]) %>%
    as.data.frame()
  
  benchmark.blood <- out %>%
    filter(case=='control',tissue=='blood') %>%
    select(disease,site,benchmean=mean,benchsd=sd,benchn=n) 
  disease.blood <- out %>%
    filter(case=='case',tissue=='blood') %>%
    select(disease,site,bloodmean=mean,bloodsd=sd,bloodn=n) 
  disease.tissue <- out %>%
    filter(case=='case',tissue=='tissue') %>%
    select(disease,site,tissuemean=mean,tissuesd=sd,tissuen=n) 
  
  out2 <- benchmark.blood %>%
    merge(disease.blood) %>%
    merge(disease.tissue) %>%
    mutate(gapblood=(bloodmean-benchmean)/benchsd) %>%
    mutate(gaptissue=(tissuemean-bloodmean)/bloodsd) %>%
    mutate(pblood=pt(gapblood,bloodn-1)) %>%
    mutate(ptissue=pt(gaptissue,tissuen-1))
  
  outi <- out2 %>%
    filter((pblood<0.025)|(pblood>0.975)) %>%
    merge(ager)
  
  results[[length(results)+1]] <- outi

  Xi <- X[rownames(X) %in% outi$site,,drop=F]
  
  datas[[length(datas)+1]] <- Xi

  chunk_start <- chunk_start + chunk_size
}

datas <- do.call(rbind,datas)
results <- do.call(rbind,results)

save(datas,results)

##########################################

results %>%
  filter(padj<0.05,abs(r)>=0.1,gapblood*gaptissue>0) %>%
  group_by(disease) %>%
  summarise(site=n_distinct(site)) %>%
  arrange(desc(site))

results %>%
  filter(padj<0.05,abs(r)>=0.1,gapblood*gaptissue>0)
