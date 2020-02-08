

setwd("~/lipids_mvp_old/data/")
i=0
file.name=paste0("BMAbatch",i,".rds")
f=readRDS(file.name)
lfsrtab=f$lfsr

for(i in 1:14){
  file.name=paste0("BMAbatch",i)
  f=readRDS(paste0(file.name,".rds"))
  lfsrtab=rbind(lfsrtab,f$lfsr)
  remove(f)
}

r=unique(rownames(lfsrtab))
c=lfsrtab[r,]
saveRDS(c,"combinedlfsrMVPBMA.rds")

###for the ukbb batches
i=0
#file.name=paste0("mashbatchUKBB",i,".rds")
#file.name=paste0("bma_ukbb_batch_V",i,".rds")
file.name=paste0("mash_ukbb_batch_V",i,".rds")
f=readRDS(file.name)
lfsrtab=f$lfsr

for(i in 1:6){
  #file.name=paste0("mashbatchUKBB",i,".rds")
  #file.name=paste0("bma_ukbb_batch_V",i,".rds")
  file.name=paste0("mash_ukbb_batch_V",i,".rds")
  f=readRDS(paste0(file.name))
  lfsrtab=rbind(lfsrtab,f$lfsr)
  remove(f)
}

r=unique(rownames(lfsrtab))
c=lfsrtab[r,]
#saveRDS(c,"combinedlfsrUKBB_bma_v.rds")
saveRDS(c,"combinedlfsrUKBB_mash_v.rds")


###for the mvp batches
i=0
file.name=paste0("bma_mvp_batch",i,".rds")
f=readRDS(file.name)
lfsrtab=f$lfsr

for(i in 1:14){
  #file.name=paste0("mashbatchUKBB",i,".rds")
  #file.name=paste0("bma_ukbb_batch_V",i,".rds")
  #file.name=paste0("mash_ukbb_batch_V",i,".rds")
  file.name=paste0("bma_mvp_batch",i,".rds")
  f=readRDS(paste0(file.name))
  lfsrtab=rbind(lfsrtab,f$lfsr)
  remove(f)
}

r=unique(rownames(lfsrtab))
c=lfsrtab[r,]
#saveRDS(c,"combinedlfsrUKBB_bma_v.rds")
#saveRDS(c,"combinedlfsrUKBB_mash_v.rds")
saveRDS(c,"combinedlfsrMVP_bma_v.rds")
rm(c);rm(lfsrtab)


i=0
#file.name=paste0("bma_mvp_batch",i,".rds")
#file.name=paste0("mash_mvp_batch",i,".rds")
#file.name=paste0("newfit_mash_ukbb_batch_V",i,".rds")
file.name=paste0("newfit_bma_ukbb_batch_V",i,".rds")
f=readRDS(file.name)
lfsrtab=f$lfsr

for(i in 1:14){
  #file.name=paste0("mashbatchUKBB",i,".rds")
  #file.name=paste0("bma_ukbb_batch_V",i,".rds")
  #file.name=paste0("mash_ukbb_batch_V",i,".rds")
  #file.name=paste0("bma_mvp_batch",i,".rds")
  #file.name=paste0("mash_mvp_batch",i,".rds")
  #file.name=paste0("newfit_mash_ukbb_batch_V",i,".rds")
  file.name=paste0("newfit_bma_ukbb_batch_V",i,".rds")
  f=readRDS(paste0(file.name))
  lfsrtab=rbind(lfsrtab,f$lfsr)
  remove(f)
}

r=unique(rownames(lfsrtab))
c=lfsrtab[r,]
#saveRDS(c,"combinedlfsrUKBB_bma_v.rds")
#saveRDS(c,"combinedlfsrUKBB_mash_v.rds")
#saveRDS(c,"combinedlfsrMVP_bma_v.rds")
#saveRDS(c,"combinedlfsrMVP_mash_v.rds")
#saveRDS(c,"combinedlfsr_ukbb_newfit_mash_v.rds")
saveRDS(c,"combinedlfsr_ukbb_newfit_bma_v.rds")
rm(c)


###To parse UKBB

c=readRDS("~/lipids_mvp/data/combinedlfsrukbb_BMA.rds")
library("tidyverse")
r=rownames(c)
t=unlist(lapply(r,function(x){
a=str_split_fixed(x,":",3)
  paste0(a[1],":",a[2])
}))

##load in the mvp bma dataset
#ukbbbma=readRDS("combinedlfsrUKBB_bma_v.rds")
#mvpbma=readRDS("~/MVPdata/combinedlfsrMVP_bma_v.rds")

ukbbbma=readRDS("~/ukbb/combinedlfsr_ukbb_newfit_bma_v.rds")
mvpbma=readRDS("~/MVPdata/combinedlfsrMVP_bma_v.rds")

#u=readRDS("~/sharednames_ukbbmvp.rds")


#names=readRDS("~/ukbb_rownames.rds")
rownames(ukbbbma)=names


thresholds=c(100,500,1000,2000,5000,10000,100000)
shared=matrix(NA,ncol=ncol(ukbbbma),nrow=length(thresholds))
colnames(shared)=colnames(ukbbbma)
rownames(shared)=thresholds
for(i in 1:ncol(ukbbbma)){
  mvp_share=mvpbma[u,i]
  ukbb_share=ukbbbma[u,i]
  o=data.frame(sort(mvp_share,decreasing = F))
  o[o<0]=0
  o2=data.frame(sort(ukbb_share,decreasing = F))
  o2[o2<0]=0
  for(j in seq(1:length(thresholds))){
    t=thresholds[j]
    a=length(intersect(rownames(o)[1:t],rownames(o2)[1:t]))
    shared[j,i]=a
  }
}

sharedbma=shared
saveRDS(sharedbma,"sharedBMA_newfit.rds")
rm(shared)

###for mash
ukbb_mash=readRDS("~/ukbb/combinedlfsr_ukbb_newfit_mash_v.rds")
mvp_mash=readRDS("~/MVPdata/combinedlfsrMVP_mash_v.rds")

names=readRDS("~/ukbb_rownames.rds")
rownames(ukbb_mash)=names

u=readRDS("~/sharednames_ukbbmvp.rds")

thresholds=c(100,500,1000,2000,5000,10000,100000)
shared=matrix(NA,ncol=ncol(ukbb_mash),nrow=length(thresholds))
colnames(shared)=colnames(ukbb_mash)
rownames(shared)=thresholds
for(i in 1:ncol(ukbb_mash)){
  mvp_share=mvp_mash[u,i]
  ukbb_share=ukbb_mash[u,i]
  o=data.frame(sort(mvp_share,decreasing = F))
  o[o<0]=0
  # o2=data.frame(sort(ukbb_share,decreasing = F))
  # o2[o2<0]=0
for(j in seq(1:length(thresholds))){
t=thresholds[j]
  #   #a=length(intersect(rownames(o)[1:t],rownames(o2)[1:t]))
a=mean(ukbb_share[rownames(o)[1:t]])###RETURN THE MEAN LFSR OF THE TOP T MVP SNPS IN UKBB
    shared[j,i]=a
  }
}


saveRDS(shared,"sharedmash_newfit.rds")
rm(shared)


##univariate approach

zukbb=as.matrix(readRDS("~/lipids_mvp_old/data/z_ukbb.rds"))
zmash=readRDS("~/lipids_mvp_old/data/zmash.rds")

zukbb=abs(zukbb)
zmash=abs(zmash)

thresholds=c(100,500,1000,2000,5000,10000,100000)
shared3=matrix(NA,ncol=ncol(zmash),nrow=length(thresholds))
colnames(shared3)=colnames(zmash)

for(i in 1:ncol(zmash)){
  ukbb_share=zukbb[u,i]
  mvp_share=zmash[u,i]
  o=data.frame(sort(mvp_share,decreasing = T))
  o[o<0]=0
  o2=data.frame(sort(ukbb_share,decreasing = T))
  o2[o2<0]=0
  for(j in seq(1:length(thresholds))){
    t=thresholds[j]
    a=length(intersect(rownames(o)[1:t],rownames(o2)[1:t]))
    shared3[j,i]=a
  }
}

saveRDS(shared3,"sharedmash_Zraw.rds")
