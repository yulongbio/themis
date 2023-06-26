auc.df=read.csv('../Data/random_split_auc.csv')
auc.df$Cohort=factor(auc.df$Cohort,levels = c('Training','Test'))
library(ggplot2)
library(reshape2)
auc_plot.df=melt(auc.df,var.ids=c('MFR','FSI','CAFF','FEM','THEMIS'),measure.vars = c(2:6),variable.name = 'Feature',value.name = 'AUC')
png('Melinus_100_AUC.png',width = 4000,height = 3000,res = 600)
ggplot(auc_plot.df,aes(x=Feature,y=AUC,fill=Cohort))+geom_boxplot()+theme_bw()+coord_cartesian(ylim = c(0, 1))
dev.off()




