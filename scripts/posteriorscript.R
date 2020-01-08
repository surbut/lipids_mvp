
library("data.table")




#zmash=readRDS("~/lipids_mvp/data/zmash.rds")
zmash=readRDS("~/ukbb/z_ukbb.rds")
#m=readRDS("~/lipids_mvp/data/mfitMVP.rds")
#m=readRDS('~/lipids_mvp/data/m_eQTL_bma.rds')
m=readRDS('~/MVPdata/mfitMVP.rds')
m=readRDS('~/MVPdata//m_eQTL_bma.rds')
#Vhat=readRDS("~/lipids_mvp/data/MVPVhat.rds")
Vhat=readRDS("~/MVPVhat.rds")

for(i in 0:13){
  start=i*2e6+1
  stop=(i+1)*2e6
  print(c(start,stop))
  library("mashr")
  mash.data=mash_set_data(zmash[start:stop,],V = Vhat,alpha = 1)
  p=mash_compute_posterior_matrices(m$fitted_g, mash.data, algorithm.version = "Rcpp")
  saveRDS(p,file = paste0("~/MVPdata/mash_mvp_batch",i,".rds"))
  #saveRDS(p,file = paste0("~/lipids_mvp/data/BMAbatch",i,".rds"))
  #saveRDS(p,file = paste0("~/ukbb/mashbatchUKBB",i,".rds"))}
}
  
i=14
start=i*2e6+1
stop=nrow(zmash)
print(c(start,stop))
library("mashr")
mash.data=mash_set_data(zmash[start:stop,],V = Vhat,alpha = 1)
p=mash_compute_posterior_matrices(m$fitted_g, mash.data, algorithm.version = "Rcpp")
saveRDS(p,file = paste0("~/lipids_mvp/data/BMAbatch",i,".rds"))

library("ashr")

for(i in 1:13){
  start=i*2e6+1
  stop=(i+1)*2e6
  print(c(start,stop))
  bhat=zmash[start:stop,1]
  sehat=rep(1,length(bhat))
  t=ashr::ash(betahat=bhat,sebetahat=sehat,g=a[[1]]$fitted_g,fixg=T)
  int=t$result$lfsr
  comb_lfsr=c(comb_lfsr,int)
}

  
  saveRDS(p,file = paste0("~/lipids_mvp/data/batch",i,".rds"))}