d = readRDS('~/simuNBMM_lmm_lmer_diff_all.rds')

plmm = readRDS('simuNBMM_plmm_log.rds')
pnbn = readRDS('simuNBMM_pnbn_log.rds')

tlmm = readRDS('./simuNBMM_tlmm.rds')
tnbn = readRDS('./simuNBMM_tnbn.rds')


differences=c(d$`Variance components` , d$Coefficients, d$`t-values`, d$`p-values`)
labels=c(rep('Variance\ncomponents', length(d$`Variance components`)), 
  rep('Coefficients',length(d$Coefficients)),
  rep('t-values', length(d$`t-values`)),
  rep('p-values', length(d$`p-values`)))
diff_df = data.frame(diff=differences, labels=labels)

ggplot(diff_df, aes(y=diff, x=labels)) +
  geom_boxplot() +
  xlab('') +
  ylab('Differences between LMM and lmer') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15, colour="black"),
        axis.text.y = element_text(size = 15, colour="black"),
        axis.title.y = element_text(size = 16, colour="black"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))



rtime = readRDS('./rtime_table_mainfig_2b.rds')
rtime_df = reshape2::melt(rtime)
colnames(rtime_df) = c('model', 'num_cell', 'time')
ggplot(rtime_df, aes(x=num_cell, y=time, group=model)) +
  geom_line(aes(linetype=model, color=model),size=0.7)+
  geom_point(aes(color=model))+scale_color_brewer(palette="Dark2")+theme_classic()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  xlab('Number of cells') +
  ylab('Computation time (min)')+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

