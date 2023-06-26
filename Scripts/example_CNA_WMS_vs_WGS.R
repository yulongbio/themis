### Plot WGS vs WMS CIN at bin level for individual example

tiff('WMS_vs_WGS_sample_CIN.tiff',width = 5000,height = 2500,res = 600)
par(mfrow=c(2,1),mar=c(4,4,1,3))
bin_1.df=read.table('../Data/example_WGS_CNA_bin.txt',sep='\t',header = T)
bin_1.df$chrom=sub('chr','',bin_1.df$chrom)
bin_1.df$chrom=as.numeric(bin_1.df$chrom)
bin_1.df=bin_1.df[order(bin_1.df$chrom,bin_1.df$start,decreasing = F),]
bin_2.df=read.table('../Data/example_WMS_CNA_bin.txt',sep='\t',header = T)
bin_2.df$chrom=sub('chr','',bin_2.df$chrom)
bin_2.df$chrom=as.numeric(bin_2.df$chrom)
bin_2.df=bin_2.df[order(bin_2.df$chrom,bin_2.df$start,decreasing = F),]
bin_1.v=paste0(bin_1.df$chrom,':',bin_1.df$start,';',bin_1.df$end)
bin_2.v=paste0(bin_2.df$chrom,':',bin_2.df$start,';',bin_2.df$end)
bin.v=union(bin_1.v,bin_2.v)
bin_chr.v=bin.v;bin_start.v=bin.v
for (i in 1:length(bin.v)){
  bin_chr.v[i]=strsplit(bin.v[i],':')[[1]][1]
  bin_start.v[i]=strsplit(strsplit(bin.v[i],':')[[1]][2],';')[[1]][1]
}
bin_all.df=data.frame(chrom=bin_chr.v,start=bin_start.v,l2rr_wgs=NA,l2rr_wms=NA)
bin_all.df$chrom=as.numeric(bin_all.df$chrom)
bin_all.df=bin_all.df[order(bin_all.df$chrom,bin_all.df$start,decreasing = F),]

sample_names=c('WMS','WGS')
bin.df=read.table('../Data/example_WMS_CNA_bin.txt',sep='\t',header = T)
bin.df$chrom=sub('chr','',bin.df$chrom)
bin.df$chrom=as.numeric(bin.df$chrom)
bin.df=bin.df[order(bin.df$chrom,bin.df$start,decreasing = F),]
baseline.m=as.matrix(bin.df[,-c(1:7)])
baseline.v=apply(baseline.m,1,mean)
log2rr=log2(bin.df$sample_norm/baseline.v)
idx=match(paste0(bin.df$chrom,':',bin.df$start),paste0(bin_all.df$chrom,':',bin_all.df$start))


bin_all.df$l2rr_wgs[idx]=log2rr
chr.v=table(bin_all.df$chrom)
gaps.v=rep(NA,length(chr.v)-1)
for (c in 2:(length(chr.v))){
  gaps.v[c-1]=sum(chr.v[1:(c-1)])
}

plot(bin_all.df$l2rr_wgs,ylim=c(-1.5,1.5),pch='.',xaxt='n',ylab='Log2 ratio',yaxt='n',frame.plot = F,col='blue',xlab='',main='WGS',cex=2)
points(which(bin_all.df$l2rr_wgs>0.3),bin_all.df$l2rr_wgs[which(bin_all.df$l2rr_wgs>0.3)],col='red',pch='.',cex=2)
points(which(bin_all.df$l2rr_wgs<(-0.3)),bin_all.df$l2rr_wgs[which(bin_all.df$l2rr_wgs<(-0.3))],col='green',pch='.',cex=2)
abline(v=c(0,gaps.v),lty='dashed')
axis(1,at=c(gaps.v,sum(chr.v))-chr.v/2,labels = names(chr.v),tick = F)
axis(2,at=c(-1.5,0,1.5),labels = c(-1.5,0,1.5))



bin.df=read.table('../Data/example_WGS_CNA_bin.txt',sep='\t',header = T)
bin.df$chrom=sub('chr','',bin.df$chrom)
bin.df$chrom=as.numeric(bin.df$chrom)
bin.df=bin.df[order(bin.df$chrom,bin.df$start,decreasing = F),]
baseline.m=as.matrix(bin.df[,-c(1:7)])
baseline.v=apply(baseline.m,1,mean)
log2rr=log2(bin.df$sample_norm/baseline.v)
idx=match(paste0(bin.df$chrom,':',bin.df$start),paste0(bin_all.df$chrom,':',bin_all.df$start))


bin_all.df$l2rr_wms[idx]=log2rr
chr.v=table(bin_all.df$chrom)
gaps.v=rep(NA,length(chr.v)-1)
for (c in 2:(length(chr.v))){
  gaps.v[c-1]=sum(chr.v[1:(c-1)])
}

plot(bin_all.df$l2rr_wms,ylim=c(-1.5,1.5),pch='.',xaxt='n',ylab='Log2 ratio',yaxt='n',frame.plot = F,col='blue',xlab='',main='WMS',cex=2)
points(which(bin_all.df$l2rr_wms>0.3),bin_all.df$l2rr_wms[which(bin_all.df$l2rr_wms>0.3)],col='red',pch='.',cex=2)
points(which(bin_all.df$l2rr_wms<(-0.3)),bin_all.df$l2rr_wms[which(bin_all.df$l2rr_wms<(-0.3))],col='green',pch='.',cex=2)
abline(v=c(0,gaps.v),lty='dashed')
axis(1,at=c(gaps.v,sum(chr.v))-chr.v/2,labels = names(chr.v),tick = F)
axis(2,at=c(-1.5,0,1.5),labels = c(-1.5,0,1.5))
dev.off()
