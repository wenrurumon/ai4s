
rm(list=ls())
setwd('/data4/ewas')
library(data.table)
library(dplyr)
library(impute)
library(reticulate)
use_python('/opt/anaconda3/bin/python')
skllm <- import('sklearn.linear_model')
np <- import('numpy')

################################################

testdata <- fread('disease/testdata.csv')
traindata <- fread('disease/traindata.csv')
raw <- fread('Rugao.Meth.beta.csv')
map <- fread('2019Rugao.Meth.Func.csv')
map2 <- fread('Sample.info.csv')
map$gender <- map2$sex

map <- map %>%
  select(id,age,mCD2019,hCD2019,SCF2019c1,gender) %>%
  melt(id=c('id','age','gender')) %>%
  group_by(id,age,gender) %>%
  summarise(value=mean(value,na.rm=T)) %>%
  # mutate(flag=ifelse(value>=0.5,1,0)) %>%
  mutate(flag=ifelse(value==1,1,0)) %>%
  filter(!is.na(age)) %>%
  select(id,age,gender,flag)
X <- as.matrix(raw[,-1])
rownames(X) <- raw$V1
raw <- t(X)
raw2 <- apply(raw,1,as.numeric)
rownames(raw2) <- colnames(raw)
raw <- log((raw2+0.0001)/(1-raw2+0.0001))
raw2 <- data.table(cpgsite=rownames(raw),raw)
ewas <- data.table(traindata,testdata[,-1])

rgsite <- raw2$cpgsite
esite <- ewas$cpgsite

raw2 <- raw2[match(esite,rgsite),]
raw2$cpgsite <- esite

################################################

ewasmap <- rbind(fread('disease/trainmap.csv'),fread('disease/testmap.csv'))[,-1] %>%
  filter(disease %in% c("Alzheimer's disease",'control')) %>%
  select(sample_id,age,gender,disease) %>%
  mutate(disease=ifelse(disease=='control',0,1))
rgmap <- map %>%
  select(sample_id=id,age,gender,disease=flag) %>%
  mutate(gender=ifelse(gender==1,'M','F'))
map <- rbind(rgmap,ewasmap)
sites <- raw2$cpgsite

############

rgage <- rgmap$age
map4sel <- lapply(0:1,function(i){
  ewasmap %>% filter(age<70,disease==i)
})
ewid <- c(map4sel[[1]]$sample_id,map4sel[[2]]$sample_id)
# ewid <- c(sample(map4sel[[1]]$sample_id,150),sample(map4sel[[2]]$sample_id,50))
ewage <- (map %>% filter(sample_id%in%ewid))$age
map <- rbind(rgmap,ewasmap %>% filter(sample_id%in%ewid)) %>%
  filter(!is.na(disease),!is.na(gender))

rgdata <- raw2 %>% select(colnames(raw2)[(colnames(raw2)%in%map$sample_id)])
ewasdata <- ewas %>% select(colnames(ewas)[colnames(ewas)%in%map$sample_id])
out <- cbind(rgdata,ewasdata) %>% as.matrix
outmap <- map[match(colnames(out),map$sample_id),]
rownames(out) <- sites

outmap <- outmap[match(sample(outmap$sample_id),outmap$sample_id),]
out <- out[,match(outmap$sample_id,colnames(out))]
outmap <- outmap %>%
  mutate(source=ifelse(substr(sample_id,1,1)=='G','rugao','ewas'))
# outmap$sample_id <- paste0('valsample',1:nrow(outmap))
# colnames(out) <- outmap$sample_id

#Alignment

test <- apply(out,2,function(x){
  c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))
})
test <- t(test)
test <- data.table(test,age=outmap$age,id=1:nrow(test)) %>%
  as.data.frame

map2align <- np$load('stat_data.npy')
map2align <- map2align %>%
  as.data.frame %>%
  select(mean=1,sd=2,age=3) %>%
  mutate(valid=1:nrow(map2align)) %>%
  as.data.frame()

test[,1:3] <- apply(test[,1:3],2,round,digits=5)
map2align[,1:3] <- apply(map2align[,1:3],2,round,digits=5)

test2 <- merge(map2align,test,all.x=T,by='age') %>%
  mutate(meangap=abs(mean.x/mean.y-1),sdgap=abs(sd.x/sd.y-1),gap=meangap+sdgap) %>%
  arrange(valid,gap)
test2 <- test2[match(unique(test2$valid),test2$valid),]

#reprocess

out2 <- out[,test2$id]
outmap2 <- outmap[test2$id,]

test <- cbind(t(apply(out2,2,function(x){c(mean(x,na.rm=T),sd(x,na.rm=T))})),outmap2$age)
outmap2$sample_id <- paste0('valsample',1:nrow(outmap2))
colnames(out2) <- outmap2$sample_id

out <- out2
outmap <- outmap2

#####################

valdata <- out
map <- outmap
valdata2 <- valdata
x <- valdata2[!is.na(valdata2)]
set.seed(123); valdata2[!is.na(valdata2)] <- 1
valmap <- map
dummymap <- valmap
dummymap$age <- 50
dummymap$disease <- 1
dummymap$gender <- 'M'

write.csv(valdata,'processeddata/valdata.csv')
write.csv(valdata2,'processeddata/dummydata.csv')
write.csv(valmap,'processeddata/valmap.csv',row.names=F)
write.csv(dummymap,'processeddata/dummymap.csv',row.names=F)
