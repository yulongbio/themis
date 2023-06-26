auc.df=read.csv('../Data/split_by_hospital_auc.csv',header = T)
library(ggplot2)
library(reshape2)
auc_train.df=auc.df[,-grep('test',colnames(auc.df))]
auc_train.df=data.frame(auc_train.df,Cohort='Training')
colnames(auc_train.df)[9:14]=c('Early_percentage','MFR','FSI','CAFF','FEM','THEMIS')
auc_test.df=auc.df[,-grep('train',colnames(auc.df))]
auc_test.df=data.frame(auc_test.df,Cohort='Test')
colnames(auc_test.df)[9:14]=c('Early_percentage','MFR','FSI','CAFF','FEM','THEMIS')
auc_split.df=rbind(auc_train.df,auc_test.df)
auc_plot.df=melt(auc_split.df,var.ids=c('MFR','FSI','CAFF','FEM','THEMIS'),measure.vars = c('MFR','FSI','CAFF','FEM','THEMIS'),variable.name = 'Feature',value.name = 'AUC')
auc_plot.df$Cohort=factor(auc_plot.df$Cohort,levels=c('Training','Test'))
# Boxplots of AUCs
png('Random_split_AUC.png',width = 3000,height = 2000,res = 600)
ggplot(auc_plot.df,aes(x=Feature,y=AUC,fill=Cohort))+geom_boxplot()+theme_bw()+coord_cartesian(ylim = c(0, 1))
dev.off()

# Plot AUC vs early_percentage
png('Random_split_prediction_vs_early_percent.png',width = 6000,height = 2500,res = 600)
ggplot(data=auc_plot.df,aes(x=Early_percentage,y=AUC))+geom_point()+theme_bw()+geom_smooth(method=lm)+facet_grid(Cohort~Feature)+xlab('Fraction of early-stage samples')+ylab('AUC')+theme(axis.title = element_text(face='bold'),strip.text = element_text(face='bold'))
dev.off()
