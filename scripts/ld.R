ukbb_ld=readRDS("~/Dropbox/combined_ld_block_ukbb.rds")
mvp=readRDS("~/Dropbox/combined_ld_block_mvp.rds")

a=which(rowSums(mvp[,c(3:6)]!=0)!=0)

b=which(rowSums(ukbb_ld[,c(3:6)]!=0)!=0)

length(a)
length(b)
length(intersect(a,b))


head(which(rowSums(ukbb_ld[,c(3:6)]!=0)==0&rowSums(mvp[,c(3:6)]!=0)!=0))
ukbb_ld[2,]
mvp[2,]

mashz=readRDS("~/Dropbox/zmash.rds")
mashz["1:851344",]

ukbb=readRDS("~/Dropbox/z_ukbb.rds")

library("tidyverse")
r=rownames(ukbb)
t=unlist(lapply(r,function(x){
  a=str_split_fixed(x,":",3)
  paste0(a[1],":",a[2])
}))
ukbb["1:1060355:G:C",]
