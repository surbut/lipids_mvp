
share=readRDS("~/Dropbox/combined_shared_lfsr_mvp_mash_chr_bp.rds")
#share=readRDS("~/Dropbox/combined_shared_lfsr_ukbb_mash_chr_bp.rds")
blocksize=500000
disco=read.csv("~/Dropbox/all_disco.csv")
d=(apply(disco,2,function(x){as.character(x)}))
thresh=0.05
s=str_split_fixed(d[,1],":",2)
##just for hdl (d[,1])
derek=data.frame(chr=paste0("chr",s[,1]),bp=as.numeric(s[,2]))
###

s=str_split_fixed(df2[,2],":",2)
derek=data.frame(chr=paste0("chr",s[,1]),bp=as.numeric(s[,2]))



for(j in 1:4){
  #s=str_split_fixed(d[,j],":",2)
  ##just for hdl (d[,1])
  #derek=data.frame(chr=paste0("chr",s[,1]),bp=as.numeric(s[,2]))
  
  
for(i in 0:21){
  mat <- data.frame(NULL,NULL,NULL)
  #lfsr_mat <- data.frame(NULL,NULL,NULL)
  chr=paste0("chr",i+1)
  df=share[share$chr==chr,]
  df$bp=as.numeric(as.character(df$b))
  start=min(df$bp) ### chromosomal boundaries
  stop=max(df$bp) 
  blocks=trunc((stop-start)/blocksize,2)
  in_chrom=df[df$chr==chr,]
  ordered_chrom=in_chrom[order(in_chrom$bp),] 
  
  derek_in_chrom=derek[derek$chr==chr,]
  ordered_derek=derek_in_chrom[order(derek_in_chrom$bp),] 
  ## return all elements of shared matrix in chromosome
  print(blocks)
  for(q in 0:blocks){
    print(c(i,q))
    pos_start=start+q*(blocksize+1)##starting position of the LD block
    pos_stop=min(pos_start+blocksize,stop)
    goodguys=ordered_chrom[ordered_chrom$bp>(pos_start-1)&ordered_chrom$bp<(pos_stop+1),]
    good_derek=ordered_derek[ordered_derek$bp>(pos_start-1)&ordered_derek$bp<(pos_stop+1),]
    row_ind=q+1
    mat[row_ind,1]=paste0(chr,":",pos_start)
    mat[row_ind,2]=paste0(chr,":",pos_stop)
    if(nrow(goodguys)==0){mat[row_ind,3]=0}
    if(nrow(goodguys)>0){
      #thresh_sat=goodguys[,c(3:6)]<0.05
      #for(j in 1:4){
        #if(colSums(thresh_sat)[j]>0){
        #mat[row_ind,3]=min(goodguys[,3])<thresh
        mat[row_ind,3]=dim(good_derek)[1]>0
      #}
    }
  }
  
  #saveRDS(mat,paste0("~/Dropbox/",chr,"hdl_min_per_block_mvp_TF.rds"))
  #saveRDS(mat,paste0("~/Dropbox/",chr,colnames(share)[j+2],"min_per_block_derek.rds"))
  saveRDS(mat,paste0("~/Dropbox/",chr,"min_per_block_cad.rds"))
}
  
  }


df=readRDS("~/Dropbox/chr1min_per_block_cad.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"min_per_block_cad.rds"));df=rbind(df,file)}




share=readRDS("~/Dropbox/combined_shared_lfsr_ukbb_mash_chr_bp.rds")
blocksize=500000
for(i in 0:21){
  mat <- data.frame(NULL,NULL,NULL)
  #lfsr_mat <- data.frame(NULL,NULL,NULL)
  chr=paste0("chr",i+1)
  df=share[share$chr==chr,]
  df$bp=as.numeric(as.character(df$b))
  start=min(df$bp) ### chromosomal boundaries
  stop=max(df$bp) 
  blocks=trunc((stop-start)/blocksize,2)
  in_chrom=df[df$chr==chr,]
  ordered_chrom=in_chrom[order(in_chrom$bp),] 
  
  # derek_in_chrom=derek[derek$chr==chr,]
  # ordered_derek=derek_in_chrom[order(derek_in_chrom$bp),] 
  # ## return all elements of shared matrix in chromosome
  print(blocks)
  for(q in 0:blocks){
    print(c(i,q))
    pos_start=start+q*(blocksize+1)##starting position of the LD block
    pos_stop=min(pos_start+blocksize,stop)
    goodguys=ordered_chrom[ordered_chrom$bp>(pos_start-1)&ordered_chrom$bp<(pos_stop+1),]
    #good_derek=ordered_derek[ordered_derek$bp>(pos_start-1)&ordered_derek$bp<(pos_stop+1),]
    row_ind=q+1
    mat[row_ind,1]=paste0(chr,":",pos_start)
    mat[row_ind,2]=paste0(chr,":",pos_stop)
    if(nrow(goodguys)==0){mat[row_ind,c(3:6)]=rep(0,4)}
    if(nrow(goodguys)>0){
    #thresh_sat=goodguys[,c(3:6)]<0.05
    for(j in 1:4){
      #if(colSums(thresh_sat)[j]>0){
      mat[row_ind,2+j]=min(goodguys[,j+2])<thresh
      #mat[row_ind,4]=dim(good_derek)[1]>0
      }
    }
  }
  #saveRDS(mat,paste0("~/Dropbox/",chr,"hdl_min_per_block_mvp_TF.rds"))
  #saveRDS(mat,paste0("~/Dropbox/",chr,"hdl_min_per_block_ukbb_TF.rds"))
  #saveRDS(mat,paste0("~/Dropbox/",chr,"500kb.rds"))
  saveRDS(mat,paste0("~/Dropbox/",chr,"500kb_ukbb.rds"))
}

f=readRDS("~/Dropbox/chr1500kb.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"500kb.rds"));f=rbind(f,file)}


f2=readRDS("~/Dropbox/chr1500kb_ukbb.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"500kb_ukbb.rds"));f2=rbind(f2,file)}

g=readRDS("~/Dropbox/encore_data/chr1500kb_encore_short.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/encore_data/chr",i,"500kb_encore_short.rds"));g=rbind(g,file)}


g2=readRDS("~/Dropbox/encore_data/chr1500kb_encore.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/encore_data/chr",i,"500kb_encore.rds"));g=rbind(g2,file)}

# 
f3=readRDS("~/Dropbox/chr1hdlmin_per_block_derek.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"hdlmin_per_block_derek.rds"));f3=rbind(f3,file)}

f4=readRDS("~/Dropbox/chr1ldlmin_per_block_derek.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"ldlmin_per_block_derek.rds"));f4=rbind(f4,file)}

f5=readRDS("~/Dropbox/chr1tgmin_per_block_derek.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"tgmin_per_block_derek.rds"));f5=rbind(f5,file)}

f6=readRDS("~/Dropbox/chr1tcmin_per_block_derek.rds")
for(i in 2:22){file=readRDS(paste0("~/Dropbox/chr",i,"tcmin_per_block_derek.rds"));f6=rbind(f6,file)}

# 
# 
# 
# 
df=data.frame(start=f3$V1,stop=f3$V2,hdl=f3$V3,ldl=f4$V3,tg=f5$V3,tc=f6$V3)
# 

venn.diagram(list(mvp=which(f$V3==1),ukbb=which(f2$V3==1),cad=which(df[,3]==1)), "~/Dropbox//cad_hdl.tiff",main="Intersection with HDL")
venn.diagram(list(mvp=which(f$V4==1),ukbb=which(f2$V4==1),cad=which(df[,3]==1)), "~/Dropbox//cad_ldl.tiff",main = "Intersection with LDL")
venn.diagram(list(mvp=which(f$V5==1),ukbb=which(f2$V5==1),cad=which(df[,3]==1)), "~/Dropbox//cad_tg.tiff",main = "Intersection with TG")
venn.diagram(list(mvp=which(f$V6==1),ukbb=which(f2$V6==1),cad=which(df[,3]==1)), "~/Dropbox//cad_tc.tiff",main = "Intersection with TC")


hdl=venn.diagram(list(mvp=which(f$V3==1),ukbb=which(f2$V3==1),previous=which(df$hdl==1)), main="Intersection with HDL","~/Dropbox//hdl.png")
ldl=venn.diagram(list(mvp=which(f$V4==1),ukbb=which(f2$V4==1),previous=which(df$ldl==1)), main = "Intersection with LDL","~/Dropbox/ldl.png")
tg=venn.diagram(list(mvp=which(f$V5==1),ukbb=which(f2$V5==1),previous=which(df$tg==1)), main = "Intersection with TG","~/Dropbox//tg.png")
tc=venn.diagram(list(mvp=which(f$V6==1),ukbb=which(f2$V6==1),previous=which(df$tc==1)), main = "Intersection with TC","~/Dropbox/tc.png")



hdl=venn.diagram(list(mvp=which(f$V3==1),ukbb=which(f2$V3==1),encore=which(g$hdl==1)), main="Intersection with HDL","~/Dropbox//hdl_es.tiff")
ldl=venn.diagram(list(mvp=which(f$V4==1),ukbb=which(f2$V4==1),encore=which(g$ldl==1)), main = "Intersection with LDL","~/Dropbox/ldl_es.png")
tg=venn.diagram(list(mvp=which(f$V5==1),ukbb=which(f2$V5==1),encore=which(g$tg==1)), main = "Intersection with TG","~/Dropbox//tg_es.png")
tc=venn.diagram(list(mvp=which(f$V6==1),ukbb=which(f2$V6==1),encore=which(g$tc==1)), main = "Intersection with TC","~/Dropbox/tc_es.png")
