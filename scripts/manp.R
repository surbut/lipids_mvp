#ukbb_mash=readRDS("~/Dropbox/combinedlfsr_ukbb_newfit_mash_v.rds")
mvp_mash=readRDS("~/Dropbox/combinedlfsrMVP_mash_v.rds")
share=readRDS("~/Dropbox/sharednames_ukbbmvp.rds")
df=readRDS("~/Dropbox/sharedchrbp.rds")
mvp=mvp_mash[share,]
#ukbb=cbind(df,ukbb_mash[share,])
#mvp=cbind(df,mvp_mash[share,])
#maxp=readRDS("~/Dropbox/max_p_mvp_raw.rds");maxp=cbind(df,maxp)

#min_ukbb=apply(ukbb[,c(3:6)],1,function(x){min(x)})
#min_mvp=apply(mvp[,c(3:6)],1,function(x){min(x)})
#mvp=cbind(df,min_mvp)
#ukbb=cbind(df,min_ukbb)
library('stringr')
numextract <- function(string){
str_extract(string, "\\-*\\d+\\.*\\d*")
}

c=numextract(df$chr)
bp=as.numeric(as.character(df$bp))
# mvp$chr=c
# mvp$bp=as.numeric(as.character(df$bp))
# mvp$min_mvp[min_mvp<0]=0

# mvp=cbind(c,as.numeric(as.character(df$bp),min_mvp)
# qqman::manhattan(mvp,chr="chr",bp="bp",p="min_mvp",genomewideline = -log10(5e-1),ylim=c(0,100))
# 
# png("mymanhattan.png", width=950, height=500)
# print(manhattan.plot(chr, pos, pvalue, sig.level=5e-8))
# dev.off()


pval=c(min_mvp,maxp$psha)
bp=c(bp,bp)
chr=c(c,c)
groups=rep( c("MASH", "naiive"), each=length(mvp$bp) )
png("mymanhattan.png", width=950, height=500)
manhattan_plot( pval, bp, chr, groups, main="MASHvsNaiive")
dev.off()


zmash=readRDS("~/Dropbox/zmash.rds")
zshare=zmash[share,]
p=apply(zshare,1,function(x){2*(1-pnorm(abs(x)))})

pval=c(mvp$hdl,p$hdl)
bp=c(bp,bp)
chr=c(c,c)
png("~/Dropbox/mymanhattan_hdl.png", width=950, height=500)
manhattan_plot( pval, bp, chr, groups, main="MASHvsNaiive_hdl")
dev.off()


pval=c(mvp$ldl,p$ldl)
bp=c(bp,bp)
chr=c(c,c)
png("~/Dropbox/mymanhattan_ldl.png", width=950, height=500)
manhattan_plot( pval, bp, chr, groups, main="MASHvsNaiive_ldl")
dev.off()


pval=c(mvp$tg,p$tg)
bp=c(bp,bp)
chr=c(c,c)
png("~/Dropbox/mymanhattan_tg.png", width=950, height=500)
manhattan_plot( pval, bp, chr, groups, main="MASHvsNaiive_tg")
dev.off()

pval=c(mvp$tc,p$tc)
bp=c(bp,bp)
chr=c(c,c)
png("~/Dropbox/mymanhattan_tc.png", width=950, height=500)
manhattan_plot( pval, bp, chr, groups, main="MASHvsNaiive_tc")
dev.off()
