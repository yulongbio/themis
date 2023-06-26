auc.df=read.csv('../Data/pc_selection_auc.csv',header = T)
library(ggplot2)
library(reshape2)
auc_train.df=auc.df[,-grep('validation',colnames(auc.df))]
auc_train.df=data.frame(auc_train.df,Cohort='Training')
colnames(auc_train.df)[2:4]=c('MFR','FSI','FEM')
auc_test.df=auc.df[,-grep('train',colnames(auc.df))]
auc_test.df=data.frame(auc_test.df,Cohort='Validation')
colnames(auc_test.df)[2:4]=c('MFR','FSI','FEM')
auc_split.df=rbind(auc_train.df,auc_test.df)
auc_plot.df=melt(auc_split.df,var.ids=c('MFR','FSI','FEM'),measure.vars = c('MFR','FSI','FEM'),variable.name = 'Feature',value.name = 'AUC')
auc_plot.df$Variance=auc_plot.df$Variance*100
auc_plot.df$Cohort=factor(auc_plot.df$Cohort,levels=c('Training','Validation'))
auc_plot.df$Variance=factor(auc_plot.df$Variance)
png('AUC_by_variance.png',width = 7500,height = 2000,res = 600)
ggplot(auc_plot.df,aes(x=Variance,y=AUC,color=Cohort))+geom_boxplot()+facet_grid(.~Feature)+xlab('Variance explained (%)')+theme_bw()
dev.off()
