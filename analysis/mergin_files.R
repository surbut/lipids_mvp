ldl100k=read.csv("ldl_Qced_unrel_100k.tsv.bgz",head=TRUE,sep='\t')
ldl=data.frame(ldl100k$locus,ldl100k$alleles,ldl100k$t_stat,ldl100k$standard_error)
hdl100k=read.csv("hdl_Qced_unrel_100k.tsv.bgz",head=TRUE,sep='\t')
hdl=data.frame(hdl100k$locus,hdl100k$alleles,hdl100k$t_stat,hdl100k$standard_error)
tg100k=read.csv("TG_LOG_Qced_unrel_100k.tsv.bgz",head=TRUE,sep='\t')
tg=data.frame(tg100k$locus,tg100k$alleles,tg100k$t_stat,tg100k$standard_error)
tc100k=read.csv("chol_Qced_unrel_100k.tsv.bgz",head=TRUE,sep='\t')
tc=data.frame(tc100k$locus,tc100k$alleles,tc100k$t_stat,tc100k$standard_error)
tstat=data.frame(hdl$hdl100k.locus,hdl$hdl100k.alleles,hdl$hdl100k.t_stat,ldl$ldl100k.t_stat,tg$tg100k.t_stat,tc$tc100k.t_stat)
sestat=data.frame(hdl$hdl100k.locus,hdl$hdl100k.alleles,hdl$hdl100k.standard_error,ldl$ldl100k.standard_error,tg$tg100k.standard_error,tc$tc100k.standard_error)
write.table(tstat,"100k_tstat_ukbb.txt",sep="\t",quote=FALSE)
##########


ldl10k=read.csv("ldl_Qced_unrel_10k.tsv.bgz",head=TRUE,sep='\t')
ldl=data.frame(ldl10k$locus,ldl10k$alleles,ldl10k$t_stat,ldl10k$standard_error)
hdl10k=read.csv("hdl_Qced_unrel_10k.tsv.bgz",head=TRUE,sep='\t')
hdl=data.frame(hdl10k$locus,hdl10k$alleles,hdl10k$t_stat,hdl10k$standard_error)
tg10k=read.csv("TG_LOG_Qced_unrel_10k.tsv.bgz",head=TRUE,sep='\t')
tg=data.frame(tg10k$locus,tg10k$alleles,tg10k$t_stat,tg10k$standard_error)
tc10k=read.csv("chol_Qced_unrel_10k.tsv.bgz",head=TRUE,sep='\t')
tc=data.frame(tc10k$locus,tc10k$alleles,tc10k$t_stat,tc10k$standard_error)
tstat=data.frame(hdl$hdl10k.locus,hdl$hdl10k.alleles,hdl$hdl10k.t_stat,ldl$ldl10k.t_stat,tg$tg10k.t_stat,tc$tc10k.t_stat)
sestat=data.frame(hdl$hdl10k.locus,hdl$hdl10k.alleles,hdl$hdl10k.standard_error,ldl$ldl10k.standard_error,tg$tg10k.standard_error,tc$tc10k.standard_error)
write.table(tstat,"10k_tstat_ukbb.txt",sep="\t",quote=FALSE)

############


ldl250k=read.csv("ldl_Qced_unrel_250k.tsv.bgz",head=TRUE,sep='\t')
ldl=data.frame(ldl250k$locus,ldl250k$alleles,ldl250k$t_stat,ldl250k$standard_error)
hdl250k=read.csv("hdl_Qced_unrel_250k.tsv.bgz",head=TRUE,sep='\t')
hdl=data.frame(hdl250k$locus,hdl250k$alleles,hdl250k$t_stat,hdl250k$standard_error)
tg250k=read.csv("TG_LOG_Qced_unrel_250k.tsv.bgz",head=TRUE,sep='\t')
tg=data.frame(tg250k$locus,tg250k$alleles,tg250k$t_stat,tg250k$standard_error)
tc250k=read.csv("chol_Qced_unrel_250k.tsv.bgz",head=TRUE,sep='\t')
tc=data.frame(tc250k$locus,tc250k$alleles,tc250k$t_stat,tc250k$standard_error)
tstat=data.frame(hdl$hdl250k.locus,hdl$hdl250k.alleles,hdl$hdl250k.t_stat,ldl$ldl250k.t_stat,tg$tg250k.t_stat,tc$tc250k.t_stat)
sestat=data.frame(hdl$hdl250k.locus,hdl$hdl250k.alleles,hdl$hdl250k.standard_error,ldl$ldl250k.standard_error,tg$tg250k.standard_error,tc$tc250k.standard_error)
write.table(tstat,"250k_tstat_ukbb.txt",sep="\t",quote=FALSE)


###########3
############


ldl50k=read.csv("ldl_Qced_unrel_50k.tsv.bgz",head=TRUE,sep='\t')
ldl=data.frame(ldl50k$locus,ldl50k$alleles,ldl50k$t_stat,ldl50k$standard_error)
hdl50k=read.csv("hdl_Qced_unrel_50k.tsv.bgz",head=TRUE,sep='\t')
hdl=data.frame(hdl50k$locus,hdl50k$alleles,hdl50k$t_stat,hdl50k$standard_error)
tg50k=read.csv("TG_LOG_Qced_unrel_50k.tsv.bgz",head=TRUE,sep='\t')
tg=data.frame(tg50k$locus,tg50k$alleles,tg50k$t_stat,tg50k$standard_error)
tc50k=read.csv("chol_Qced_unrel_50k.tsv.bgz",head=TRUE,sep='\t')
tc=data.frame(tc50k$locus,tc50k$alleles,tc50k$t_stat,tc50k$standard_error)
tstat=data.frame(hdl$hdl50k.locus,hdl$hdl50k.alleles,hdl$hdl50k.t_stat,ldl$ldl50k.t_stat,tg$tg50k.t_stat,tc$tc50k.t_stat)
sestat=data.frame(hdl$hdl50k.locus,hdl$hdl50k.alleles,hdl$hdl50k.standard_error,ldl$ldl50k.standard_error,tg$tg50k.standard_error,tc$tc50k.standard_error)
write.table(tstat,"50k_tstat_ukbb.txt",sep="\t",quote=FALSE)

#################
