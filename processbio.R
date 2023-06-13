
rm(list=ls())
library(data.table)
library(dplyr)
setwd('/data4/ewas/disease')
load('disease_methylation_v1.RData')
load('sample_disease.RData')

rawmap <- sample_disease
rawmap2 <- t(disease_download[1:4,])
raw <- t(disease_download[-1:-4,])
raw2 <- apply(raw,1,as.numeric)
rownames(raw2) <- colnames(raw)
raw <- log((raw2+0.0001)/(1-raw2+0.0001))

map <- rawmap %>% select(sample_id,age,gender=sex,sample_type,disease)
map %>% 
    group_by(disease,sample_type) %>%
    summarise(n=n(),age=mean(age,na.rm=T)) %>%
    as.data.frame()

pdata <- data.table(map,t(raw))
pdata <- pdata %>% filter(!is.na(age))
pdata$disease <- ifelse(is.na(pdata$disease),'control',pdata$disease)

test <- pdata %>%
    filter(!is.na(age)) %>%
    group_by(disease,sample_type) %>%
    summarise(n=n(),mean=mean(age,na.rm=T),sd=sd(age,na.rm=T)) %>%
    mutate(z=(mean-48.8)/25.5) %>%
    arrange(desc(z)) %>%
    filter(z>0) 
pdata <- pdata %>% filter(!is.na(age),disease%in%test$disease) 
plist <- lapply(unique(pdata$disease),function(d){
    pdata %>% filter(disease==d)
})
names(plist) <- unique(pdata$disease)

train_test <- lapply(plist,function(x){
    sel <- sample(1:nrow(x))/nrow(x)
    list(train=x[(sel)<0.8,,drop=F],test=x[(sel)>0.8,,drop=F])
})
trainset <- do.call(rbind,lapply(train_test,function(x){x$train}))
testset <- do.call(rbind,lapply(train_test,function(x){x$test}))
trainset$sample_id <- paste0('train',1:nrow(trainset)+10000)
testset$sample_id <- paste0('test',1:nrow(testset)+10000)

trainmap <- trainset[,1:5]
write.csv(trainmap,'trainmap.csv')
testmap <- testset[,1:5]
write.csv(testmap,'testmap.csv')

testdata <- t(testset[,-1:-5])
colnames(testdata) <- testmap$sample_id
fwrite(testdata,'testdata.csv',sep=',')

traindata <- t(trainset[,-1:-5])
colnames(traindata) <- trainmap$sample_id
fwrite(traindata,'traindata.csv',sep=',')

############################################################
############################################################


