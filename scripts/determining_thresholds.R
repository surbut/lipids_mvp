pmat=readRDS("pmat.rds")
ash=readRDS("ash_zmash.rds")

hdl=cbind(sort(pmat[,"hdl"]),sort(ash[,1]))
ldl=cbind(sort(pmat[,"ldl"]),sort(ash[,2]))
tg=cbind(sort(pmat[,"tg"]),sort(ash[,3]))
hdl[310,]
ldl[840,]
tg[10,]