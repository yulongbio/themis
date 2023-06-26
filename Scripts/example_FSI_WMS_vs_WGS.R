tiff('WMS_vs_WGS_sample_FSI.tiff',width = 5000,height = 2500,res = 600)

par(mfrow=c(2,1),mar=c(4,4,1,3))
data.df=read.table('../Data/example_WGS_WMS_FSI_bin.txt',header = T)
data_z.df=data.df
for (i in 1:ncol(data.df)){
  data_z.df[,i]=(data.df[,i]-mean(data.df[,i]))/sd(data.df[,i])
}
regions.df=read.table('../Data/FSI_regions.txt',header = T)
chr.v=sub('chr','',regions.df$seqnames)
chr.v=as.numeric(chr.v)
start.v=regions.df$start
chr_start.df=data.frame(chr=chr.v,start=start.v)
idx=with(chr_start.df,order(chr,start,decreasing = F))
data_z.df=data_z.df[idx,]
chr_name.v=table(chr.v)
gaps.v=rep(NA,length(chr_name.v)-1)
for (c in 2:(length(chr_name.v))){
  gaps.v[c-1]=sum(chr_name.v[1:(c-1)])
}
healthy_wms.m=readRDS('../Data/healthy_baseline_FSI.rds')[['WGS']]
plot(data_z.df[,1],type='l',main='WGS',ylim=c(-2,4),xaxt='n',ylab='Z score',yaxt='n',frame.plot = F,col='black',xlab='')
for (i in 1:ncol(healthy_wgs.m)){
  lines(healthy_wgs.m[,i],col='grey60')
}
lines(data_z.df[,1],lwd=2,col='black')

abline(v=c(0,gaps.v),lty='dashed')
axis(1,at=c(gaps.v,sum(chr_name.v))-chr_name.v/2,labels = names(chr_name.v),tick = F)
axis(2,at=c(-2,0,4),labels = c(-2,0,4))
plot(data_z.df[,2],type = 'l',main='WMS',ylim=c(-2,4),xaxt='n',ylab='Z score',yaxt='n',frame.plot = F,col='black',xlab='')
healthy_wgs.m=readRDS('../Data/healthy_baseline_FSI.rds')[['WMS']]

for (i in 1:ncol(healthy_wms.m)){
  lines(healthy_wms.m[,i],col='grey60')
}
lines(data_z.df[,2],lwd=2,col='black')
abline(v=c(0,gaps.v),lty='dashed')
axis(1,at=c(gaps.v,sum(chr_name.v))-chr_name.v/2,labels = names(chr_name.v),tick = F)
axis(2,at=c(-2,0,4),labels = c(-2,0,4))
dev.off()
