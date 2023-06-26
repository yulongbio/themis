library(tidyverse)
library(stats)
library(openxlsx)

ci_plot <- function(file_in){
    file_in <- file_in
    file_out <- 'THEMIS_sensitivity.hist.pdf'
    df <- read.table(file_in, header = T)
    df$Cohort<-sub('Discovery','Training',df$Cohort)
    df$Cohort<-sub('Validation','Test',df$Cohort)
    df$Cohort<-factor(df$Cohort,levels = c('Training','Test'))
    df=df[order(df$Cohort,decreasing = F),]
    df <- df[df$Group != 'OV',]
    df <- df[df$Group != 'HEALTHY',]
    df <- df[df$Stage != 'unknown',]
    df <- df[!is.na(df$Stage),]
    df$Group <- as.character(df$Group)
    df$Group <- ifelse(df$Group == 'Cancer', 'OVERALL', df$Group)
    df$CI_low <- as.numeric(gsub("-.*", "", df$CI))
    df$CI_high <- as.numeric(gsub(".*-", "", df$CI))
    df$Sens <- df$Sens
    df$CI_low <- df$CI_low
    df$CI_high <- df$CI_high
    df <- df %>% group_by(Stage, Group) %>% mutate(group_number = paste0(total, collapse=' | '))
    df$Stage_new <- paste0(df$Stage, "\n", '(', df$group_number ,')')
    df <- transform(df, Group = factor(Group, levels = c('BRCA', 'COREAD', 'ESCA', 'LIHC', 'NSCLC', 'PACA','STAD','OVERALL')))
    Stage_new <- df[df$Stage_new == 'Training',]
    pdf(file_out, width = 8, height = 9,pointsize = 11)
    p <- ggplot(data = df, aes(x = Stage_new, y = Sens, ymin = CI_low, ymax = CI_high, colour = Cohort, fill = Cohort)) +
         geom_histogram(stat = "identity", width = 0.6,  position = position_dodge(0.7)) + 
         geom_errorbar(position = position_dodge(width = 0.7), width = 0.2, color = 'black',size=0.4) + 
             scale_color_manual(values = c("#27B7B7", "#FF66B2")) + 
             scale_fill_manual(values = c("#27B7B7", "#FF66B2")) +
         geom_point(position = position_dodge(width = 0.7), size = 0.5, color = 'black') + 
         theme_bw() +
         facet_wrap(vars(Group), ncol = 3, scales = "free_x") +
         theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),legend.position = 'none') +
         scale_y_continuous(limits  = c(0, 100), labels = function(x) paste0(x, "%"), breaks = seq(0, 100, by = 25)) + 
         ylab('Sensitivity')  + 
         xlab('Stage')
    print(p)
    dev.off()
}

ci_plot('../Data/THEMIS_sensitivity.txt')

