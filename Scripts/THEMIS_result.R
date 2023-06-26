data.df=read.table('../Data/MONITOR_prediction_scores.txt',header = T,check.names = F,na.strings = 'AAA')
data.df=data.df[,-c(2,3)]
colnames(data.df)=c('Patient_ID','Diagnosis','Stage','Cohort','MFR','FSI','CAFF','FEM','THEMIS')
data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

validation.df=subset(data.df,Cohort=='Test')
discovery.df=subset(data.df,Cohort=='Training')

cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')

### Pairwise correlation between omics
library(PerformanceAnalytics)
chart.Correlation<-function (R, histogram = TRUE, method = c("pearson", "kendall", "spearman"), ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method = cormeth, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram) 
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel,...)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,...)
}

png('omics_pairwise_cor.png',width = 4000,height = 4000,res = 600,pointsize = 12)
chart.Correlation(data.df[,c('MFR','FSI','CAFF','FEM')],histogram = TRUE,pch='.',method = 'spearman')
dev.off()

########Plot THEMIS ROC by cancer 
plot_roc_by_cancer<-function(op_num,comp,cohort){
  if(cohort=='Training'){data.df=discovery.df}
  if(cohort=='Test'){data.df=validation.df}
  op_names=c('FSI','CAFF','MFR','FEM','THEMIS')
  data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))
  cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')
  png(paste0(op_names[op_num],'_',cohort,'_by_cancer_ROC.png'),width = 3950,height = 3950,res = 600,pointsize = 14)
  library(pROC)
  library(RColorBrewer)
  plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[2]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.4-0.05*1,print.auc.x=0.5,main=cohort)
  text(0.51,0.4-0.05*1,labels = paste0(cancer_type.v[2],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[1])
  for (i in 3:length(cancer_type.v)){
    plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[i]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[i-1],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.4-0.05*(i-1),print.auc.x=0.5)
    text(0.51,0.4-0.05*(i-1),labels = paste0(cancer_type.v[i],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[i-1])
  }
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  dev.off()
}
plot_roc_by_cancer(5,'<','Training')
plot_roc_by_cancer(5,'<','Test')

###### Plot THEMIS ROC by stage
plot_roc_by_phase<-function(op_num,comp,cohort){
  if(cohort=='Training'){data.df=discovery.df}
  if(cohort=='Test'){data.df=validation.df}
  op_names=c('FSI','CAFF','MFR','FEM','THEMIS')
  data.df$Stage=factor(data.df$Stage,levels = c('HEALTHY','I','II','III','IV','NA'))
  Phase.v=c('HEALTHY','I','II','III','IV','NA')
  png(paste0(op_names[op_num],'_',cohort,'_by_stage_ROC.png'),width = 3950,height = 3950,res = 600,pointsize = 14)
  library(pROC)
  library(RColorBrewer)
  plot.roc(data.df$Stage,data.df[,op_names[op_num]],levels=c('HEALTHY',Phase.v[2]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(Phase.v),'Dark2')[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.35-0.05*1,print.auc.x=0.5,main=cohort)
  text(0.52,0.35-0.05*1,labels = paste0(Phase.v[2],' '),adj = c(1,1),col=brewer.pal(length(Phase.v),'Dark2')[1])
  for (i in 3:length(Phase.v)){
    plot.roc(data.df$Stage,data.df[,op_names[op_num]],levels=c('HEALTHY',Phase.v[i]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(Phase.v),'Dark2')[i-1],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.35-0.05*(i-1),print.auc.x=0.5)
    text(0.52,0.35-0.05*(i-1),labels = paste0(Phase.v[i],''),adj = c(1,1),col=brewer.pal(length(Phase.v),'Dark2')[i-1])
  }
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  plot.roc(data.df$cancer,data.df[,op_names[op_num]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col='black',type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.35-0.05*(i),print.auc.x=0.5)
  text(0.52,0.35-0.05*(i),labels = paste0('OVERALL'),adj = c(1,1),col='black')
  dev.off()
}
plot_roc_by_phase(5,'<','Training')
plot_roc_by_phase(5,'<','Test')


######### Plot ROC for individual modalities
plot_roc_by_cancer_train_test<-function(op_num,comp){
  library(pROC)
  library(RColorBrewer)
  op_names=c('FSI','CAFF','MFR','FEM')
  tiff(paste0(op_names[op_num],'_by_cancer_ROC_train_test.tif'),width = 6000,height = 3000,res = 600,pointsize = 10.5,compression = 'lzw')
  par(mfrow=c(1,2),mar=c(4,4,4,4))
  data.df=discovery.df
  data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))
  cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')
  
  plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[2]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.45-0.05*1,print.auc.x=0.45,main='Training')
  text(0.452,0.45-0.05*1,labels = paste0(cancer_type.v[2],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[1])
  for (i in 3:length(cancer_type.v)){
    plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[i]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[i-1],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.45-0.05*(i-1),print.auc.x=0.45)
    text(0.452,0.45-0.05*(i-1),labels = paste0(cancer_type.v[i],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[i-1])
  }
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  plot.roc(data.df$cancer,data.df[,op_names[op_num]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col='black',type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.45-0.05*i,print.auc.x=0.45)
  text(0.452,0.45-0.05*i,labels = paste0('OVERALL '),adj = c(1,1),col='black')
  
  
  data.df=validation.df
  data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))
  cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')
  
  plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[2]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.45-0.05*1,print.auc.x=0.45,main='Test')
  text(0.452,0.45-0.05*1,labels = paste0(cancer_type.v[2],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[1])
  for (i in 3:length(cancer_type.v)){
    plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[i]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[i-1],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.45-0.05*(i-1),print.auc.x=0.45)
    text(0.452,0.45-0.05*(i-1),labels = paste0(cancer_type.v[i],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[i-1])
  }
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  plot.roc(data.df$cancer,data.df[,op_names[op_num]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col='black',type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.45-0.05*i,print.auc.x=0.45)
  text(0.452,0.45-0.05*i,labels = paste0('OVERALL '),adj = c(1,1),col='black')
  dev.off()
}


for (i in 1:4){
  plot_roc_by_cancer_train_test(i,'<')
}

### Boxplots of prediction scores by cancer
library(RColorBrewer)
library(ggplot2)
library(patchwork)
cutoff.v=rep(NA,5)
for (i in 1:5){
  h_train.df=subset(discovery.df,Diagnosis=='HEALTHY')
  idx=round(nrow(h_train.df)*0.99)
  cutoff.v[i]=sort(h_train.df[,4+i],decreasing = F)[idx]
}
names(cutoff.v)=colnames(discovery.df[5:9])

op_names=c('FSI','MFR','FEM','THEMIS')
for (i in op_names){
png(paste0(i,'_scores_by_cancer.png'),width = 4000,height = 4000,res = 600)
p1=ggplot(discovery.df,aes_string(x='Diagnosis',y=i,color='Diagnosis'))+geom_boxplot(outlier.shape = NA,position=position_dodge(width = 0.3),width=0.6)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(7,'Dark2')))+ylim(c(0,1))+labs(x="",y='Predictive score',title = 'Training')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')
p2=ggplot(validation.df,aes_string(x='Diagnosis',y=i,color='Diagnosis'))+geom_boxplot(outlier.shape = NA,position=position_dodge(width = 0.3),width=0.6)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(7,'Dark2')))+ylim(c(0,1))+labs(x="Diagnosis",y='Predictive score',title = 'Test')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'bottom')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')+scale_fill_discrete(name = "Diagnosis")
plot(p1/p2)
dev.off()
}

i='CAFF'
png(paste0(i,'_scores_by_cancer.png'),width = 4000,height = 4000,res = 600)
p1=ggplot(discovery.df,aes_string(x='Diagnosis',y=i,color='Diagnosis'))+geom_boxplot(outlier.shape = NA,position=position_dodge(width = 0.3),width=0.6)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(7,'Dark2')))+ylim(c(0,35))+labs(x="",y='PA score',title = 'Training')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')
p2=ggplot(validation.df,aes_string(x='Diagnosis',y=i,color='Diagnosis'))+geom_boxplot(outlier.shape = NA,position=position_dodge(width = 0.3),width=0.6)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(7,'Dark2')))+ylim(c(0,35))+labs(x="Diagnosis",y='PA score',title = 'Test')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'bottom')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')+scale_fill_discrete(name = "Diagnosis")
plot(p1/p2)
dev.off()

### Boxplots of prediction scores by stage
library(RColorBrewer)
library(ggplot2)
library(patchwork)
cutoff.v=rep(NA,5)
for (i in 1:5){
  h_train.df=subset(discovery.df,Diagnosis=='HEALTHY')
  idx=round(nrow(h_train.df)*0.99)
  cutoff.v[i]=sort(h_train.df[,4+i],decreasing = F)[idx]
}
names(cutoff.v)=colnames(discovery.df[5:9])

op_names=c('FSI','MFR','FEM','THEMIS')
for (i in op_names){
  tiff(paste0(i,'_score_by_stage.tiff'),width = 3000,height = 4000,res = 600,compression = 'lzw')
  p1=ggplot(discovery.df,aes_string(x='Stage',y=i,color='Stage'))+geom_boxplot(outlier.shape = NA,position=position_dodge(width = 0.3),width=0.6)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(5,'Dark2')))+ylim(c(0,1))+labs(x="",y='Predictive score',title = 'Training')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')
  p2=ggplot(validation.df,aes_string(x='Stage',y=i,color='Stage'))+geom_boxplot(outlier.shape = NA,position=position_dodge(width = 0.3),width=0.6)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(5,'Dark2')))+ylim(c(0,1))+labs(x="Stage",y='Predictive score',title = 'Test')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'bottom')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')+scale_fill_discrete(name = "Stage")
  plot(p1/p2)
  dev.off()
}

i='CAFF'
tiff(paste0(i,'_scores_by_stage.tiff'),width = 3000,height = 4000,res = 600,compression = 'lzw')
p1=ggplot(discovery.df,aes_string(x='Stage',y=i,color='Stage'))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(7,'Dark2')))+ylim(c(0,35))+labs(x="",y='PA score',title = 'Training')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')
p2=ggplot(validation.df,aes_string(x='Stage',y=i,color='Stage'))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=20,position = position_jitter(0.2))+theme_classic()+scale_color_manual(values = c('black',brewer.pal(7,'Dark2')))+ylim(c(0,35))+labs(x="Stage",y='PA score',title = 'Test')+theme(plot.title = element_text(hjust = 0.5),legend.position = 'bottom')+geom_hline(yintercept = cutoff.v[i],linetype='dashed')+scale_fill_discrete(name = "Stage")
plot(p1/p2)
dev.off()



##### Plot overall ROC for all models
plot_roc_overall_combine_omics_with_ensemble<-function(comp,cohort){
  if(cohort=='Training'){data.df=discovery.df}
  if(cohort=='Test'){data.df=validation.df}
  op_names=c('MFR','FSI','CAFF','FEM','THEMIS')
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  tiff(paste0('Overall_',cohort,'_omics_with_ensemble_ROC.tif'),width = 3400,height = 3400,res = 600,compression = 'lzw',pointsize = 13)
  library(pROC)
  library(RColorBrewer)
  cols=c(brewer.pal(length(op_names),'Dark2'))
  plot.roc(data.df$cancer,data.df[,op_names[1]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col=cols[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.3-0.05*1,print.auc.x=0.6,lwd=3,main=cohort)
  text(0.61,0.3-0.05*1,labels = paste0(op_names[1],' '),adj = c(1,1),col=cols[1])
  for (i in 2:5){
    plot.roc(data.df$cancer,data.df[,op_names[i]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col=cols[i],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.3-0.05*i,print.auc.x=0.6,lwd=3)
    text(0.61,0.3-0.05*i,labels = paste0(op_names[i],' '),adj = c(1,1),col=cols[i])
    
  }
  dev.off()
}
plot_roc_overall_combine_omics_with_ensemble('<','Training')
plot_roc_overall_combine_omics_with_ensemble('<','Test')


### Calculate cancer sensitivity for THEMIS at 99% specificity
data.df=read.table('all_cohort_scores.txt',header = T,check.names = F,na.strings = 'AAA')
data.df=data.df[,-c(2,3)]
colnames(data.df)=c('PID','Diagnosis','Stage','Cohort','MFR','FSI','CAFF','FEM','THEMIS')
data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))
data.df$CAFF=10^(data.df$CAFF)
data.df=cbind(data.df,cancer='OVERALL')
data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
cutoff.v=rep(NA,5)
for (i in 1:5){
  h_train.df=subset(discovery.df,Diagnosis=='HEALTHY')
  idx=round(nrow(h_train.df)*0.99)
  cutoff.v[i]=sort(h_train.df[,4+i],decreasing = F)[idx]
}
names(cutoff.v)=colnames(discovery.df[5:9])

pos.l=list()
for (i in 1:5){
  data_cancer.df=subset(data.df,Diagnosis!='HEALTHY')
  pos.l[[i]]=data_cancer.df[which(data_cancer.df[,4+i]>=cutoff.v[i]),]$PID
}
names(pos.l)=colnames(discovery.df)[5:9]
library(RColorBrewer)
library(venn)
library(ggplot2)
library(ggpolypath)
png('Overlap_0.99_spec_cancer_training_test.png',res=600,height = 3500,width = 3500,pointsize = 12)
venn(list(`MFR`=pos.l[['MFR']],`FSI`=pos.l[['FSI']],`CAFF`=pos.l[['CAFF']],`FEM`=pos.l[['FEM']],`THEMIS`=pos.l[['THEMIS']]),zcolor = brewer.pal(5,'Dark2')[c(1,2,3,4,5)],sncs = 1.4,ilcs=1.1,box = F)
dev.off()


#Calculate healthy specificity for THEMIS at 99% specificity
pos.l=list()
for (i in 1:5){
  data_cancer.df=subset(data.df,Diagnosis=='HEALTHY')
  pos.l[[i]]=data_cancer.df[which(data_cancer.df[,4+i]>=cutoff.v[i]),]$PID
}
names(pos.l)=colnames(discovery.df)[5:9]

library(RColorBrewer)
library(venn)
library(ggplot2)
library(ggpolypath)
png('Overlap_0.99_spec_healthy_training_test.png',res=600,height = 3500,width = 3500,pointsize = 12)
venn(list(`MFR`=pos.l[['MFR']],`FSI`=pos.l[['FSI']],`CAFF`=pos.l[['CAFF']],`FEM`=pos.l[['FEM']],`THEMIS`=pos.l[['THEMIS']]),zcolor = brewer.pal(5,'Dark2')[c(1,2,3,4,5)],sncs = 1.4,ilcs=1.1,box = F)
dev.off()







