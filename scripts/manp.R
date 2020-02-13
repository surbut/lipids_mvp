#ukbb_mash=readRDS("~/Dropbox/combinedlfsr_ukbb_newfit_mash_v.rds")
mvp_mash=readRDS("~/Dropbox/combinedlfsrMVP_mash_v.rds")
share=readRDS("~/Dropbox/sharednames_ukbbmvp.rds")
df=readRDS("~/Dropbox/sharedchrbp.rds")

ukbb=cbind(df,ukbb_mash[share,])
mvp=cbind(df,mvp_mash[share,])
maxp=readRDS("~/Dropbox/max_p_mvp_raw.rds");maxp=cbind(df,maxp)

min_ukbb=apply(ukbb[,c(3:6)],1,function(x){min(x)})
min_mvp=apply(mvp[,c(3:6)],1,function(x){min(x)})
mvp=cbind(df,min_mvp)
ukbb=cbind(df,min_ukbb)
numextract <- function(string){
str_extract(string, "\\-*\\d+\\.*\\d*")
}
c=numextract(df$chr)
mvp$chr=c
mvp$bp=as.numeric(as.character(df$bp))
mvp$min_mvp[min_mvp<0]=0

mvp=cbind(c,as.numeric(as.character(df$bp),min_mvp)
qqman::manhattan(mvp,chr="chr",bp="bp",p="min_mvp",genomewideline = -log10(5e-1),ylim=c(0,100))

png("mymanhattan.png", width=950, height=500)
print(manhattan.plot(chr, pos, pvalue, sig.level=5e-8))
dev.off()


pval=c(min_mvp,maxp$psha)
bp=c(mvp$bp,mvp$bp)
chr=c(mvp$chr,mvp$chr)
groups=rep( c("MASH", "naiive"), each=length(mvp$bp) )
png("mymanhattan.png", width=950, height=500)
manhattan_plot( pval, bp, chr, groups, main="MASHvsNaiive")
dev.off()