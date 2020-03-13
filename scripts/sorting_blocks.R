
share=readRDS("~/Dropbox/combined_shared_lfsr_mvp_mash_chr_bp.rds")
mat <- data.frame(NULL,NULL,NULL)
for(i in 0:21){
  mat <- data.frame(NULL,NULL,NULL)
  chr=paste0("chr",i+1)
  df=share[share$chr==chr,]
  df$bp=as.numeric(as.character(df$b))
  start=min(df$bp) ### chromosomal boundaries
  stop=max(df$bp) 
  blocks=trunc((stop-start)/blocksize,2)
  in_chrom=df[df$chr==chr,]
  ordered_chrom=in_chrom[order(in_chrom$bp),] ## return all elements of shared matrix in chromosome
  print(blocks)
  for(q in 0:blocks){
    print(c(i,q))
    pos_start=start+q*(blocksize+1)##starting position of the LD block
    pos_stop=min(pos_start+blocksize,stop)
    goodguys=ordered_chrom[ordered_chrom$bp>(pos_start-1)&ordered_chrom$bp<(pos_stop+1),]
    row_ind=q+1
    mat[row_ind,1]=paste0(chr,":",pos_start)
    mat[row_ind,2]=paste0(chr,":",pos_stop)
    if(nrow(goodguys)==0){mat[row_ind,c(3:6)]=rep(0,4)}
    if(nrow(goodguys)>0){
    thresh_sat=goodguys[,c(3:6)]<thresh
    for(j in 1:4){
        if(colSums(thresh_sat)[j]>0){
          mat[row_ind,2+j]=rownames(goodguys)[which.min(goodguys[,2+j])]
        } else {
          mat[row_ind,2+j]=0
          }
      }
    }
  }
  saveRDS(mat,paste0("~/Dropbox/",chr,"min_per_block.rds"))
}

> sum(rowSums(UKBB[,c(3:6)]!=0)!=0)
[1] 7697
> sum(rowSums(UKBB[,c(3:6)]!=0)==0)
[1] 48037
> dim(UKBB)
[1] 55734     6
> MVP=data.frame(MVP)
> sum(rowSums(MVP[,c(3:6)]!=0)!=0)
[1] 24280

##########for stroing lfsr


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
  ordered_chrom=in_chrom[order(in_chrom$bp),] ## return all elements of shared matrix in chromosome
  print(blocks)
  for(q in 0:blocks){
    print(c(i,q))
    pos_start=start+q*(blocksize+1)##starting position of the LD block
    pos_stop=min(pos_start+blocksize,stop)
    goodguys=ordered_chrom[ordered_chrom$bp>(pos_start-1)&ordered_chrom$bp<(pos_stop+1),]
    row_ind=q+1
    mat[row_ind,1]=paste0(chr,":",pos_start)
    mat[row_ind,2]=paste0(chr,":",pos_stop)
    if(nrow(goodguys)==0){mat[row_ind,c(3:6)]=rep(0,4)}
    if(nrow(goodguys)>0){
      #thresh_sat=goodguys[,c(3:6)]<0.05
      for(j in 1:4){
        #if(colSums(thresh_sat)[j]>0){
        mat[row_ind,2*j+1]=rownames(goodguys)[which.min(goodguys[,2+j])]
        mat[row_ind,2*j+2]=max(min(goodguys[,2+j]),0)
      }
    }
  }
  saveRDS(mat,paste0("~/Dropbox/",chr,"min_per_block_all_loci.rds"))
  #saveRDS(lfsr_mat,paste0("~/Dropbox/",chr,"lfsr_min_per_block_all_loci.rds"))
  }


      