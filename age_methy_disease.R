
rm(list=ls())
setwd('/mnt/Shared_Data/ewas')
library(reticulate)
use_python('/home/huzixin/anaconda3/envs/test/bin/python')
library(data.table)
library(dplyr)
np <- import('numpy')
skl <- import('sklearn')

#Import data

raw <- fread("methylation_disease_wholeblood.csv")
map <- fread('map_disease_wholeblood2.csv')[,-1:-2]
X <- raw[,-1] %>% as.matrix %>% t
colnames(X) <- raw[[1]]

rm(raw)
gc()
X <- X[,colMeans(is.na(X))<0.5]

#Age Associate

x <- X[is.na(map$disease),,drop=F]
xage <- map$age[is.na(map$disease)]

X.cor <- cor(x,xage,use='pairwise')
temp <- x; temp[is.na(temp)] <- 0
temp.cor <- cor(temp,xage,use='pairwise')
X <- X[,(abs(X.cor)>0.2)&(abs(temp.cor)>0.2)]
Z <- map %>%
  select(age,gender=sex,disease) %>%
  mutate(disease=tolower(disease)) %>%
  mutate(gender=(gender=='M')+0) %>%
  mutate(disease=ifelse(is.na(disease),'control',disease))

X <- X[!is.na(Z$age),]
Z <- Z[!is.na(Z$age),]
Z %>%
  group_by(disease) %>%
  summarise(n(),min(age),max(age),mean(age),sd(age))

d2model <- Z %>%
  group_by(disease) %>%
  summarise(mean=mean(age),sd=sd(age)) %>%
  filter(disease!='control') %>%
  data.table(
    Z %>%
      filter(disease=='control') %>%
      summarise(benchmean=mean(age),benchsd=sd(age))
  ) %>%
  mutate(z=(mean-benchmean)/benchsd,pval=pnorm(z)) %>%
  mutate(pval=ifelse(pval>0.5,(1-pval)*2,pval*2)) %>%
  filter(z>1)

X <- X[Z$disease %in% c('control',d2model$disease),,drop=F]
Z <- Z[Z$disease %in% c('control',d2model$disease),]

#Benchmark

X2 <- X
X2[is.na(X2)] <- 0
X2 <- log((X2+0.0001)/(1-X2+0.0001),base=2)
X.bench <- lapply(unique(Z$disease),function(d){
  data.table(
    disease=d,
    site=colnames(X2),
    mean=colMeans(X2[Z$disease==d,],na.rm=T),
    sd=apply(X2[Z$disease==d,],2,sd,na.rm=T)
  )
})
names(X.bench) <- unique(Z$disease)
X.bench <- lapply(X.bench[names(X.bench)!='control'],function(x){
  x %>% 
    merge(X.bench$control %>% select(site,benchmean=mean,benchsd=sd)) %>%
    mutate(z=(mean-benchmean)/benchsd,pval=pnorm(z)) %>%
    mutate(pval=ifelse(pval>0.5,(1-pval)*2,pval*2))
})
X.bench <- do.call(rbind,X.bench)

X.bench <- X.bench %>%
  as.data.frame() %>%
  merge(data.table(site=row.names(X.cor),
    cor=as.numeric(X.cor),
    cor2=as.numeric(temp.cor)))

X.candidate <- X.bench %>% 
  filter(pval<0.05,cor*z>0) %>%
  group_by(disease) %>%
  summarise(n=n())

#Model

d2model <- unique((X.bench %>% filter(pval<0.05,cor*z>0))$disease)
models <- lapply(d2model,function(d){
  print(d)
  x <- X2[Z$disease%in%c(d,'control'),
    (X.bench %>% filter(disease==d,pval<0.05,cor*z>0))$site,drop=F]
  y <- ((Z$disease[Z$disease%in%c(d,'control')])!='control')+0
  model <- skl$linear_model$LogisticRegression(max_iter=10000L)
  model$fit(scale(x),y)
  coef_ <- model$coef_
  rlt <- sapply(1:49,function(q){
    q <- q/50
    model$fit(x[,abs(coef_)>=quantile(abs(coef_),q),drop=F],y)
    p <- model$predict(x[,abs(coef_)>=quantile(abs(coef_),q),drop=F])
    out1 <- sum((p==0)&(y==0))
    out2 <- sum((p==0)&(y==1))
    out3 <- sum((p==1)&(y==0))
    out4 <- sum((p==1)&(y==1))
    c(q=q,
      features=sum(abs(coef_)>=quantile(abs(coef_),q)),
      out=c(out1,out2,out3,out4))
  }) %>% t
  rlt <- data.table(disease=d,rlt) %>%
    mutate(precision=out4/sum(y==1))
  rlt
})
models <- do.call(rbind,models)
models <- models %>%
  merge(
    models %>%
      merge(
        models %>% 
          group_by(disease) %>% 
          summarise(precision=max(precision)) %>%
          filter(precision>0.8)
      ) %>%
      group_by(disease) %>%
      summarise(features=min(features))
  ) %>%
  mutate(n=out3+out4)

#Finalize Model

i <- 1
d <- models$disease[i]
q <- models$q[i]
x <- X2[Z$disease%in%c(d,'control'),
  (X.bench %>% filter(disease==d,pval<0.05,cor*z>0))$site,drop=F]
y <- ((Z$disease[Z$disease%in%c(d,'control')])!='control')+0
z <- Z$age[Z$disease%in%c(d,'control')]
model <- skl$linear_model$LogisticRegression(max_iter=10000L)
model$fit(scale(x),y)
coef_ <- model$coef_
x <- x[,abs(coef_)>=quantile(abs(coef_),q),drop=F]
model$fit(x,y)
p <- model$predict(x)
table(p,y)

mvalue <- x
bvalue <- X[rownames(x),colnames(x)]
bvalue[is.na(bvalue)] <- 0
addvalue <- t(
  sapply(1:nrow(mvalue),function(i){
    coef(lm(bvalue~z))[2,]
}))

b2m <- function(x){
  x[x<0] <- 0
  x[x>1] <- 1
  log((x+0.0001)/(1-x+0.0001),base=2)
}

test <- data.table(
  Z[Z$disease%in%c(d,'control'),],
  sapply(0:30,function(i){
    model$predict(b2m(bvalue+addvalue*i))
  }) 
)

test %>%
  filter(disease=='control') %>%
  group_by(V30) %>%
  summarise(n(),mean(age),sd(age))
  
