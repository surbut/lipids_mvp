zukbb=cbind("c"=a[,2],readRDS("~/Dropbox/z_ukbb_shared_chr_bp.rds"))
zukbb$bp=as.numeric(as.character(zukbb$bp))
zukbb=zukbb[with(zukbb,order(c,bp)),]
z2=apply(zukbb[,c(4:7)],2,function(x){as.numeric(x)})
pukbb=data.frame(2*pnorm(-1*abs(z2)))
ukbb=data.frame("chr"=zukbb$chr,"bp"=zukbb$bp,"hdl"=pukbb$hdl,"ldl"=pukbb$ldl,"tg"=pukbb$tg,"tc"=pukbb$tc)

zmvp=cbind("c"=a[,2],readRDS("~/Dropbox/z_mvp_shared_chr_bp.rds"))
zmvp$bp=as.numeric(as.character(zmvp$bp))
zmvp=zmvp[with(zmvp,order(c,bp)),]
z3=apply(zmvp[,c(4:7)],2,function(x){as.numeric(x)})
pmvp=data.frame(2*pnorm(-1*abs(z3)))
mvp=data.frame("chr"=zmvp$chr,"bp"=zmvp$bp,"hdl"=pmvp$hdl,"ldl"=pmvp$ldl,"tg"=pmvp$tg,"tc"=pmvp$tc)

boundaries=read.table("~/Dropbox/block_boundaries.txt",header=T)


thresholds=c(5e-8,5e-6,5e-4,5e-3,5e-2)
traits=c("hdl","ldl","tg","tc")

marray=array(NA,dim=c(length(thresholds),nrow(boundaries),length(traits)))
uarray=array(NA,dim=c(length(thresholds),nrow(boundaries),length(traits)))


for(t in 1:length(thresholds)){
  thresh=thresholds[t]
  
  for(i in 1:nrow(boundaries)){
    
    chrom=boundaries[i,1]
    start=boundaries[i,2]
    stop=boundaries[i,3]
    
    in_chrom_mvp=mvp[mvp$chr==chrom,]
    in_chrom_ukbb=ukbb[ukbb$chr==chrom,]
    #block.start=min(in_chrom_mvp$bp[in_chrom_mvp$bp>(start-1)])
    #block.stop=max(in_chrom_mvp$bp[in_chrom_mvp$bp<(stop+1)])
    good_mvp=in_chrom_mvp[in_chrom_mvp$bp>(start-1)&in_chrom_mvp$bp<(stop+1),]
    
    if(nrow(good_mvp)==0){
      marray[t,i,]=rep(0,4)
      uarray[t,i,]=rep(0,4)
    }
    if(nrow(good_mvp)>0){
      good_ukbb=in_chrom_ukbb[in_chrom_ukbb$bp>(start-1)&in_chrom_ukbb$bp<(stop+1),]
      
      for(r in 1:length(traits)){
        lipid=traits[r]
        marray[t,i,r]=min(good_mvp[lipid])<thresh
        uarray[t,i,r]=min(good_ukbb[lipid])<thresh
        print(c(t,i,r))
      }
    }
  }
}

uarray=readRDS("~/Dropbox/rawz_marray.rds")
marray=readRDS("~/Dropbox/rawz_uarray.rds")

thresholds=c(5e-8,5e-6,5e-4,5e-3,5e-2)

op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

for(i in 1:2){
  
  a=(uarray[i,,])
  b=(marray[i,,])
  int=sapply(seq(1:4),function(x)
  {length(intersect(which(a[,x]==1),which(b[,x]==1)))
  })
  ukb_rep=round(int/colSums(a),2)
  mvp_rep=round(int/colSums(b),2)
  
  mat=matrix(c(colSums(a),colSums(b)),byrow=T,nrow=2)
  x<-barplot(mat, beside=TRUE,horiz=FALSE,#names=traits,
             ylim=c(0,max(mat)+500),
             legend.text=c("ukbb","mvp"),args.legend=list(bty="n"),
             col=c("aquamarine3","coral"),border="white",
             xlab="Lipid Traits",main=paste("RawZ, Pval",thresholds[i]))
  text(x,mat+250,matrix(c(ukb_rep,mvp_rep),byrow = T,nr=2),las=2)
}

for(i in 3:4){
  
  a=(uarray[i,,])
  b=(marray[i,,])
  int=sapply(seq(1:4),function(x)
  {length(intersect(which(a[,x]==1),which(b[,x]==1)))
  })
  ukb_rep=round(int/colSums(a),2)
  mvp_rep=round(int/colSums(b),2)
  
  mat=matrix(c(colSums(a),colSums(b)),byrow=T,nrow=2)
  x<-barplot(mat, beside=TRUE,horiz=FALSE,names=traits,ylim=c(0,max(mat)+500),
             col=c("aquamarine3","coral"),border="white",
             xlab="Lipid Traits")
  text(x,mat+250,matrix(c(ukb_rep,mvp_rep),byrow = T,nr=2),las=2)}
}
#x<-as.matrix(citysales[,2:4])

#text(x+2,y,labels=as.character(x))



ukbb_lfsr=readRDS("~/Dropbox/ukbb_refit_lfsr_shared.rds")
mvp_lfs=readRDS("~/Dropbox/combined_shared_lfsr_mvp_mash_chr_bp.rds")

mvp=cbind("c"=a[,2],mvp_lfs)
ukb=cbind("c"=a[,2],ukbb_lfsr)

ukb$bp=as.numeric(as.character(ukb$bp))
ukb=ukb[with(ukb,order(c,bp)),]
ukbb_lf=data.frame("chr"=ukb$chr,"bp"=ukb$bp,"hdl"=ukb$hdl,"ldl"=ukb$ldl,"tg"=ukb$tg,"tc"=ukb$tc)

mvp$bp=as.numeric(as.character(mvp$bp))
mvp=mvp[with(mvp,order(c,bp)),]
mvp_lf=data.frame("chr"=mvp$chr,"bp"=mvp$bp,"hdl"=mvp$hdl,"ldl"=mvp$ldl,"tg"=mvp$tg,"tc"=mvp$tc)

mvp=mvp_lf
ukbb=ukbb_lf



uarray=readRDS("~/Dropbox/uarray_lfsr.rds")
marray=readRS("~/Dropbox/marray_lfsr.rds")
thresholds=c(5e-3,5e-2,1e-1,5e-1)

marray=array(NA,dim=c(length(thresholds),nrow(boundaries),length(traits)))
uarray=array(NA,dim=c(length(thresholds),nrow(boundaries),length(traits)))


for(t in 1:length(thresholds)){
  thresh=thresholds[t]
  
  for(i in 1:nrow(boundaries)){
    
    chrom=boundaries[i,1]
    start=boundaries[i,2]
    stop=boundaries[i,3]
    
    in_chrom_mvp=mvp[mvp$chr==chrom,]
    in_chrom_ukbb=ukbb[ukbb$chr==chrom,]
    #block.start=min(in_chrom_mvp$bp[in_chrom_mvp$bp>(start-1)])
    #block.stop=max(in_chrom_mvp$bp[in_chrom_mvp$bp<(stop+1)])
    good_mvp=in_chrom_mvp[in_chrom_mvp$bp>(start-1)&in_chrom_mvp$bp<(stop+1),]
    
    if(nrow(good_mvp)==0){
      marray[t,i,]=rep(0,4)
      uarray[t,i,]=rep(0,4)
    }
    if(nrow(good_mvp)>0){
      good_ukbb=in_chrom_ukbb[in_chrom_ukbb$bp>(start-1)&in_chrom_ukbb$bp<(stop+1),]
      
      for(r in 1:length(traits)){
        lipid=traits[r]
        marray[t,i,r]=min(good_mvp[lipid])<thresh
        uarray[t,i,r]=min(good_ukbb[lipid])<thresh
        print(c(t,i,r))
      }
    }
  }
}

saveRDS(uarray,"~/Dropbox/uarray_lfsr.rds")
saveRDS(marray,"~/Dropbox/marray_lfsr.rds")


op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.2,
          mar = c(0,0,1,1) + 0.2)

uarray=readRDS("~/Dropbox/uarray_lfsr.rds")
marray=readRDS("~/Dropbox/marray_lfsr.rds")
thresholds=c(5e-3,5e-2,1e-1,5e-1)


for(i in 1:length(thresholds)){
  a=(uarray[i,,])
  b=(marray[i,,])
  int=sapply(seq(1:4),function(x)
  {length(intersect(which(a[,x]==1),which(b[,x]==1)))
  })
  
  ukb_rep=round(int/colSums(a),2)
  mvp_rep=round(int/colSums(b),2)
  
  mat=matrix(c(colSums(a),colSums(b)),byrow=T,nrow=2)
  x<-barplot(mat, beside=TRUE,horiz=FALSE,names=traits,ylim=c(0,max(mat)+500),
             legend.text=c("ukbb","mvp"),args.legend=list(bty="n"),
             col=c("aquamarine3","coral"),border="white",
             xlab="Lipid Traits",main=paste("RawZ,LFSR",thresholds[i]))
  text(x,mat+100,matrix(c(ukb_rep,mvp_rep),byrow = T,nr=2))}
