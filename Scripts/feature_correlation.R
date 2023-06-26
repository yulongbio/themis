### Plot md vs fsi vs cnv heatmap
samples.df=read.csv('../Data/feature_correlation_samples.csv',header = T)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(preprocessCore)

cnv.df=read.table('../Data/CIN_matrix.txt',header = T,check.names = F)
cnv.m=as.matrix(cnv.df[,-1])
health.v=samples.df[which(samples.df$Group=='HEALTHY'),]$Sample
coread.v=samples.df[which(samples.df$Group=='COREAD'),]$Sample
cnv.m=cnv.m[,c(health.v,coread.v)]
na.idx=which(apply(cnv.m,1,anyNA))
cnv.m=cnv.m[-na.idx,]
cnv_z.m=cnv.m
for (i in 1:ncol(cnv_z.m)){
  t.m=as.matrix(cbind(cnv.m[,1:20],cnv.m[,i]))
  t.m=normalize.quantiles(t.m)
  for (t in 1:ncol(t.m)){
    t.m[,t]=(t.m[,t]-mean(t.m[,t]))/sd(t.m[,t])
  }
  z.v=rep(NA,nrow(t.m))
  for (n in 1:nrow(t.m)){
    a=t.m[n,21]
    b=t.m[n,1:20]
    z.v[n]=(a-mean(b))/sd(b)
  }
  cnv_z.m[,i]=z.v
}
annotation_col = data.frame(T=c(rep('HEALTHY',20),rep('COREAD',20)))
annotation_col$T=factor(annotation_col$T,levels=c('HEALTHY','COREAD'))
rownames(annotation_col) = c(health.v,coread.v)
cnv_ph<-pheatmap(cnv_z.m,cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-2,2,length.out = 7),show_rownames = F,show_colnames = F,main='CNA',annotation_col = annotation_col,legend = F,annotation_legend = F,cellwidth = 4,cellheight = 0.25)
h_cnv_high.v=rep(NA,nrow(cnv_z.m));h_cnv_low.v=rep(NA,nrow(cnv_z.m))
c_cnv_high.v=rep(NA,nrow(cnv_z.m));c_cnv_low.v=rep(NA,nrow(cnv_z.m))
for (i in 1:nrow(cnv_z.m)){
  h_cnv_high.v[i]=length(which(cnv_z.m[i,1:20]>2))/20
  h_cnv_low.v[i]=length(which(cnv_z.m[i,1:20]<(-2)))/20
  c_cnv_high.v[i]=length(which(cnv_z.m[i,21:40]>2))/20
  c_cnv_low.v[i]=length(which(cnv_z.m[i,21:40]<(-2)))/20
}


md.df=read.table('../Data/MFR_matrix.txt',header = T,check.names = F)
md.m=as.matrix(md.df[,-1])
health.v=samples.df[which(samples.df$Group=='HEALTHY'),]$PID
coread.v=samples.df[which(samples.df$Group=='COREAD'),]$PID
md.m=md.m[,c(health.v,coread.v)]
md.m=md.m[-na.idx,]
md_z.m=md.m
for (i in 1:ncol(md_z.m)){
  t.m=as.matrix(cbind(md.m[,1:20],md.m[,i]))
  t.m=normalize.quantiles(t.m)
  for (t in 1:ncol(t.m)){
    t.m[,t]=(t.m[,t]-mean(t.m[,t]))/sd(t.m[,t])
  }
  z.v=rep(NA,nrow(t.m))
  for (n in 1:nrow(t.m)){
    a=t.m[n,21]
    b=t.m[n,1:20]
    z.v[n]=(a-mean(b))/sd(b)
  }
  md_z.m[,i]=z.v
}


annotation_col = data.frame(T=c(rep('HEALTHY',20),rep('COREAD',20)))
annotation_col$T=factor(annotation_col$T,levels=c('HEALTHY','COREAD'))
rownames(annotation_col) = c(health.v,coread.v)
md_ph<-pheatmap(md_z.m,cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-2,2,length.out = 7),show_rownames = F,show_colnames = F,main='MFR',annotation_col = annotation_col,legend = F,annotation_legend = F,cellwidth = 4,cellheight = 0.25)
h_md_high.v=rep(NA,nrow(md_z.m));h_md_low.v=rep(NA,nrow(md_z.m))
c_md_high.v=rep(NA,nrow(md_z.m));c_md_low.v=rep(NA,nrow(md_z.m))

for (i in 1:nrow(md_z.m)){
  h_md_high.v[i]=length(which(md_z.m[i,1:20]>2))/20
  h_md_low.v[i]=length(which(md_z.m[i,1:20]<(-2)))/20
  c_md_high.v[i]=length(which(md_z.m[i,21:40]>2))/20
  c_md_low.v[i]=length(which(md_z.m[i,21:40]<(-2)))/20
}

fsi.df=read.table('../Data/fsi_ratio_matrix.txt',header = T)
health.v=samples.df[which(samples.df$Group=='HEALTHY'),]$Sample
coread.v=samples.df[which(samples.df$Group=='COREAD'),]$Sample
rownames(fsi.df)=fsi.df$sample_id
fsi_h.m=t(as.matrix(fsi.df[health.v,-1]))
fsi_c.m=t(as.matrix(fsi.df[coread.v,-1]))
fsi.m=cbind(fsi_h.m,fsi_c.m)
fsi.m=fsi.m[-na.idx,]
fsi_z.m=fsi.m
for (i in 1:ncol(fsi_z.m)){
  t.m=as.matrix(cbind(fsi.m[,1:20],fsi.m[,i]))
  t.m=normalize.quantiles(t.m)
  for (t in 1:ncol(t.m)){
    t.m[,t]=(t.m[,t]-mean(t.m[,t]))/sd(t.m[,t])
  }
  z.v=rep(NA,nrow(t.m))
  for (n in 1:nrow(t.m)){
    a=t.m[n,21]
    b=t.m[n,1:20]
    z.v[n]=(a-mean(b))/sd(b)
  }
  fsi_z.m[,i]=z.v
}


annotation_col = data.frame(T=c(rep('HEALTHY',20),rep('COREAD',20)))
annotation_col$T=factor(annotation_col$T,levels=c('HEALTHY','COREAD'))
rownames(annotation_col) = c(health.v,coread.v)
fsi_ph<-pheatmap(fsi_z.m,cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-2,2,length.out = 7),show_rownames = F,show_colnames = F,main='FSI',annotation_col = annotation_col,legend = F,annotation_legend = F,cellwidth = 4,cellheight = 0.25)
h_fsi_high.v=rep(NA,nrow(fsi_z.m));h_fsi_low.v=rep(NA,nrow(fsi_z.m))
c_fsi_high.v=rep(NA,nrow(fsi_z.m));c_fsi_low.v=rep(NA,nrow(fsi_z.m))

for (i in 1:nrow(fsi_z.m)){
  h_fsi_high.v[i]=length(which(fsi_z.m[i,1:20]>2))/20
  h_fsi_low.v[i]=length(which(fsi_z.m[i,1:20]<(-2)))/20
  c_fsi_high.v[i]=length(which(fsi_z.m[i,21:40]>2))/20
  c_fsi_low.v[i]=length(which(fsi_z.m[i,21:40]<(-2)))/20
}

bin.m=matrix(0,nrow(fsi_z.m),6)
colnames(bin.m)=c('h_FSI','h_MD','h_CNV','c_FSI','c_MD','c_CNV')
bin.m[which(h_fsi_high.v>0.5),1]=1;bin.m[which(h_fsi_low.v>0.5),1]=(-1)
bin.m[which(h_md_high.v>0.5),2]=1;bin.m[which(h_md_low.v>0.5),2]=(-1)
bin.m[which(h_cnv_high.v>0.5),3]=1;bin.m[which(h_cnv_low.v>0.5),3]=(-1)
bin.m[which(c_fsi_high.v>0.5),4]=1;bin.m[which(c_fsi_low.v>0.5),4]=(-1)
bin.m[which(c_md_high.v>0.5),5]=1;bin.m[which(c_md_low.v>0.5),5]=(-1)
bin.m[which(c_cnv_high.v>0.5),6]=1;bin.m[which(c_cnv_low.v>0.5),6]=(-1)
annotation_col = data.frame(T=c(rep('HEALTHY',3),rep('COREAD',3)))
annotation_col$T=factor(annotation_col$T,levels=c('HEALTHY','COREAD'))

rownames(annotation_col) = colnames(bin.m)

bin_ph<-pheatmap(bin.m,cluster_rows = F,cluster_cols = F,color = c('blue','white','red'),breaks = seq(-2,2,length.out = 4),show_rownames = F,show_colnames = F,main='Differential bins',annotation_col = annotation_col,legend = T,annotation_legend = F,cellwidth = 20,cellheight = 0.25)
tiff('bin_heatmap.tif',width = 6500,height = 5000,pointsize = 12,compression = 'lzw',res = 600)
grid.arrange(fsi_ph[[4]],md_ph[[4]],cnv_ph[[4]],bin_ph[[4]],ncol=4)
dev.off()

#plot bin legend
png('bin_legend.png',width = 3000,height = 2000,res = 600)
t=seq(0,1,length.out=20)
plot(t)
legend(2,0.8,bty='n',legend = c('z > 2','z < -2'),fill = c('red','blue'),border = 'white')
dev.off()
## Plot heatmap
tiff('heatmap_legend.tif',width = 4000,height = 5000,pointsize = 12,compression = 'lzw',res = 600)
annotation_col = data.frame(T=c(rep('HEALTHY',20),rep('COREAD',20)))
annotation_col$T=factor(annotation_col$T,levels=c('HEALTHY','COREAD'))
rownames(annotation_col) = c(health.v,coread.v)
fsi_ph<-pheatmap(fsi_z.m,cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-2,2,length.out = 7),show_rownames = F,show_colnames = F,main='FSI',annotation_col = annotation_col,legend = T,annotation_legend = T,cellwidth = 4,cellheight = 0.25)
dev.off()


tiff('fsi_vs_mfr.tif',width = 3000,height = 2000,pointsize = 12,compression = 'lzw',res = 600)
cor.v=rep(NA,20)
for (i in 1:20){
  cor.v[i]=cor(fsi_z.m[,20+i],md_z.m[,20+i],method = 'spearman')
}
hist(cor.v,breaks = seq(-0.8,0.8,length.out = 101),main = 'FSI vs MFR',xlab = '',xlim=c(-0.8,0.8))
dev.off()

tiff('fsi_vs_cnv.tif',width = 3000,height = 2000,pointsize = 12,compression = 'lzw',res = 600)
cor.v=rep(NA,20)
for (i in 1:20){
  cor.v[i]=cor(fsi_z.m[,20+i],cnv_z.m[,20+i],method = 'spearman')
}
hist(cor.v,breaks = seq(-0.8,0.8,length.out = 101),main = 'FSI vs CNA',xlab = '',xlim=c(-0.8,0.8))
dev.off()


tiff('mfr_vs_cnv.tif',width = 3000,height = 2000,pointsize = 12,compression = 'lzw',res = 600)
cor.v=rep(NA,20)
for (i in 1:20){
  cor.v[i]=cor(md_z.m[,20+i],cnv_z.m[,20+i],method = 'spearman')
  
}
hist(cor.v,breaks = seq(-0.8,0.8,length.out = 101),main = 'MFR vs CNA',xlab = '',xlim=c(-0.8,0.8))
dev.off()
