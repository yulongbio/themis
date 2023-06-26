### Plot heat maps of cluster p-values
library(pheatmap)
library(RColorBrewer)

sample.df=read.csv('../Data//positive_samples_for_CSO_100spec.csv',header = T)
sample_train.df=subset(sample.df,Cohort=='Training')


data.m=readRDS('../Data/CSO_feature.rds')
# Short fragments 
mat1=matrix(NA,18,7)
rownames(mat1)=paste0('Cluster_',1:18)
colnames(mat1)=c('BRCA','NSCLC','ESCA','STAD','COREAD','LIHC','PACA')
for (cancer_type in colnames(mat1)){
  target.v=sample_train.df$Patient_ID[which(sample_train.df$Group==cancer_type)]
  other.v=sample_train.df$Patient_ID[which(sample_train.df$Group!=cancer_type)]
  
    for (row_name in rownames(mat1)){
      mat1[row_name,cancer_type]=wilcox.test(as.numeric(data.m[row_name,target.v]),as.numeric(data.m[row_name,other.v]),alternative='less')[['p.value']]
    }
      
}
bk<-seq(0,10,by=0.1)
cluster_labels=c('01 - Kidney / Bile duct','02 - Colon','03 - Breast (non-basal)','04 - Prostate','05 - Brain','06 - Thyroid','07 - Skin', '08 - Squamous','09 - Liver','10 - Nerve cell','11 - Testicular','12 - Lung (adeno)','13 - Digestive','14 - Breast (basal)','15 - Uterine','16 - Bladder','17 - Mesothelium','18 - Adrenal')
pheatmap(-log10(mat1),color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),breaks = bk,cluster_rows = F,cluster_cols = F,labels_row = cluster_labels,angle_col = '45',cellheight = 15,cellwidth = 20,fontsize=10,filename = 'short_18_p_value.pdf',main='Short fragments',legend=F,width = 8,height = 6,show_rownames = F)
# Long fragments
mat1=matrix(NA,18,7)
rownames(mat1)=paste0('Cluster_',1:18,'_long')
colnames(mat1)=c('BRCA','NSCLC','ESCA','STAD','COREAD','LIHC','PACA')
for (cancer_type in colnames(mat1)){
  target.v=sample_train.df$Sample[which(sample_train.df$Group==cancer_type)]
  other.v=sample_train.df$Sample[which(sample_train.df$Group!=cancer_type)]
  
  for (row_name in rownames(mat1)){
    mat1[row_name,cancer_type]=wilcox.test(as.numeric(data.m[row_name,target.v]),as.numeric(data.m[row_name,other.v]),alternative='less')[['p.value']]
  }
  
}
bk<-seq(0,10,by=0.1)
cluster_labels=c('01 - Kidney / Bile duct','02 - Colon','03 - Breast (non-basal)','04 - Prostate','05 - Brain','06 - Thyroid','07 - Skin', '08 - Squamous','09 - Liver','10 - Nerve cell','11 - Testicular','12 - Lung (adeno)','13 - Digestive','14 - Breast (basal)','15 - Uterine','16 - Bladder','17 - Mesothelium','18 - Adrenal')
pheatmap(-log10(mat1),color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),breaks = bk,cluster_rows = F,cluster_cols = F,labels_row = cluster_labels,angle_col = '45',cellheight = 15,cellwidth = 20,fontsize=10,filename = 'long_18_p_value.pdf',main='Long fragments',legend=F,width = 8,height = 6,show_rownames = T)
# Methylation
mat1=matrix(NA,18,7)
rownames(mat1)=paste0('Cluster_',1:18,'_hypo')
colnames(mat1)=c('BRCA','NSCLC','ESCA','STAD','COREAD','LIHC','PACA')
for (cancer_type in colnames(mat1)){
  target.v=sample_train.df$Sample[which(sample_train.df$Group==cancer_type)]
  other.v=sample_train.df$Sample[which(sample_train.df$Group!=cancer_type)]
  
  for (row_name in rownames(mat1)){
    mat1[row_name,cancer_type]=wilcox.test(as.numeric(data.m[row_name,target.v]),as.numeric(data.m[row_name,other.v]),alternative='less')[['p.value']]
  }
  
}
bk<-seq(0,10,by=0.1)
cluster_labels=c('01 - Kidney / Bile duct','02 - Colon','03 - Breast (non-basal)','04 - Prostate','05 - Brain','06 - Thyroid','07 - Skin', '08 - Squamous','09 - Liver','10 - Nerve cell','11 - Testicular','12 - Lung (adeno)','13 - Digestive','14 - Breast (basal)','15 - Uterine','16 - Bladder','17 - Mesothelium','18 - Adrenal')
pheatmap(-log10(mat1),color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),breaks = bk,cluster_rows = F,cluster_cols = F,labels_row = cluster_labels,angle_col = '45',cellheight = 15,cellwidth = 20,fontsize=10,filename = 'methylation_18_p_value.pdf',main='Methylation',legend=F,width = 8,height = 6,show_rownames = F)



### Plot confusion matrices
library(ggplot2)
library(pheatmap)
library(reshape2)
Group_sort<-c('BRCA','NSCLC','ESCA','STAD','COREAD','LIHC','PACA')
mat<-read.delim('../Data/CSO_test_prediction.txt',header=T,stringsAsFactors = F)
Group_top1<-dcast(as.data.frame(table(mat$Group,mat$TOP1)),Var1 ~ Var2)
class(Group_top1)
Group_top1
rownames(Group_top1)<-Group_top1$Var1
Group_top1<-Group_top1[,-1]
Group_top1<-Group_top1[Group_sort,Group_sort]
value<-Group_top1
Group_top1$sum<-apply(Group_top1,1,function(x){sum(x)})
Group_top1
tmp<-apply(Group_top1,2,function(x){x/Group_top1$sum})
tmp<-tmp[,-8]

bk<-seq(0,1,by=0.25)
pheatmap(tmp,display_numbers = value ,fontsize_number = 12,cluster_rows = F,cluster_cols = F,angle_col = '45',cellheight = 30,cellwidth = 30,fontsize=10,filename = 'Training_CSO.pdf',breaks=seq(0,1,length.out = 100),legend_breaks = seq(0,1.25,0.25),legend_labels = seq(0,1.25,0.25)*100,main='Training')



Group_sort<-c('BRCA','NSCLC','ESCA','STAD','COREAD','LIHC','PACA')
mat<-read.delim('../Data/CSO_test_prediction.txt',header=T,stringsAsFactors = F)
Group_top1<-dcast(as.data.frame(table(mat$Group,mat$TOP1)),Var1 ~ Var2)
class(Group_top1)
Group_top1
rownames(Group_top1)<-Group_top1$Var1
Group_top1<-Group_top1[,-1]
Group_top1<-Group_top1[Group_sort,Group_sort]
value<-Group_top1
Group_top1$sum<-apply(Group_top1,1,function(x){sum(x)})
Group_top1
tmp<-apply(Group_top1,2,function(x){x/Group_top1$sum})
tmp<-tmp[,-8]

bk<-seq(0,1,by=0.25)
pheatmap(tmp,display_numbers = value ,fontsize_number = 12,cluster_rows = F,cluster_cols = F,angle_col = '45',cellheight = 30,cellwidth = 30,fontsize=10,filename = 'Testing_CSO.pdf',breaks=seq(0,1,length.out = 100),legend_breaks = seq(0,1.25,0.25),legend_labels = seq(0,1.25,0.25)*100,main='testing')
