### Plot AUC
library(readxl)
auc.df=read_excel('../Data/downsample_AUC.xlsx')
library(ggplot2)
library(reshape2)
data.df=melt(auc.df,var.ids=c('MFR','FSI','CAFF','FEM','THEMIS'),measure.vars = c(3:7),variable.name = 'Feature',value.name = 'AUC')
data.df$Depth=factor(data.df$Depth,levels = c('3 M','15 M','30 M','60 M','90 M','120 M'))
data.df$Cohort=factor(data.df$Cohort,levels = c('Training','Test'))
png('Depth_AUC.png',width = 5000,height = 3000,res = 600)

ggplot(data.df,aes(x=Depth,y=AUC,fill=Cohort))+geom_bar(position='dodge',stat='identity')+facet_wrap(~Feature,scales = 'free_x')+theme_bw()
dev.off()

