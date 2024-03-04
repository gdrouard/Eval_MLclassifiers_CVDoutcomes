## dfstoreind is a dataframe with individual data and predictions (probabilities to belong to a class)
## dfstore_mod includes summary statistics of ML classifiers for each omic
# X0: being <1sd from the mean. X1: being within 1sd from the mean. X2: being >1sd from the mean


### below are plots that can be generated to display model performance
library(ggplot2)

#### Probability of not deviating from +1sd of the mean in non +1sd-deviant individuals (test cohort)
dfstore_ind2 <- dfstore_ind[dfstore_ind$Trueclass != "X0" &
                              dfstore_ind$Trueclass != "X1",]

dfstore_ind2$notbeX2 <- 1-dfstore_ind2$probX2

dfstore_ind2$OmicsModeling <- "Methylation CpG signature"
dfstore_ind2[dfstore_ind2$Omics== "Transcriptomics" &
               dfstore_ind2$Supervision== "SSAE","OmicsModeling"] <- "Transcriptomics (SSAE)"
dfstore_ind2[dfstore_ind2$Omics== "Transcriptomics" &
               dfstore_ind2$Supervision== "USAE","OmicsModeling"] <- "Transcriptomics (USAE)"
dfstore_ind2[dfstore_ind2$Omics== "Metabolomics" &
               dfstore_ind2$Supervision== "SSAE","OmicsModeling"] <- "Metabolomics (SSAE)"
dfstore_ind2[dfstore_ind2$Omics== "Metabolomics" &
               dfstore_ind2$Supervision== "USAE","OmicsModeling"] <- "Metabolomics (USAE)"
dfstore_ind2[dfstore_ind2$Omics== "Multi-omics" &
               dfstore_ind2$Supervision== "None","OmicsModeling"] <- "Multi-omics"


dfstore_ind2$Model <- as.factor(dfstore_ind2$Model)
ggplot(dfstore_ind2, aes(target, notbeX2, fill= OmicsModeling)) +
  geom_boxplot(aes(group_by=Omics),outlier.shape=" ",alpha=0.9,col="black") +
  theme_classic() +
  scale_fill_manual(values=c("#481568FF" ,"#453781FF","darkgrey","#FDE725FF",
                             "#287D8EFF","#1F968BFF")) +
  ylab("Probability of not deviating from +1sd of the mean in non +1sd-deviant individuals (test cohort)") +
  scale_y_continuous(limits = c(0.7, 1)) +
  coord_flip() +
  facet_wrap(~Model)  + xlab("") +
  geom_hline(aes(yintercept=0.9))


### Probability of deviating from +1sd of the mean in +1sd-deviant individuals (test cohort)
dfstore_ind2 <- dfstore_ind[dfstore_ind$Trueclass == "X2",]

dfstore_ind2$OmicsModeling <- "Methylation CpG signature"
dfstore_ind2[dfstore_ind2$Omics== "Transcriptomics" &
               dfstore_ind2$Supervision== "SSAE","OmicsModeling"] <- "Transcriptomics (SSAE)"
dfstore_ind2[dfstore_ind2$Omics== "Transcriptomics" &
               dfstore_ind2$Supervision== "USAE","OmicsModeling"] <- "Transcriptomics (USAE)"
dfstore_ind2[dfstore_ind2$Omics== "Metabolomics" &
               dfstore_ind2$Supervision== "SSAE","OmicsModeling"] <- "Metabolomics (SSAE)"
dfstore_ind2[dfstore_ind2$Omics== "Metabolomics" &
               dfstore_ind2$Supervision== "USAE","OmicsModeling"] <- "Metabolomics (USAE)"
dfstore_ind2[dfstore_ind2$Omics== "Multi-omics" &
               dfstore_ind2$Supervision== "None","OmicsModeling"] <- "Multi-omics"


dfstore_ind2$Model <- as.factor(dfstore_ind2$Model)
ggplot(dfstore_ind2, aes(target, probX2, fill= OmicsModeling)) +
  geom_boxplot(aes(group_by=Omics),outlier.shape=" ",alpha=0.9,col="black") +
  theme_classic() +
  scale_fill_manual(values=c("#481568FF" ,"#453781FF","darkgrey","#FDE725FF",
                             "#287D8EFF","#1F968BFF")) +
  ylab("Probability of deviating from +1sd of the mean in +1sd-deviant individuals (test cohort)") +
  scale_y_continuous(limits = c(0, 1)) +
  coord_flip() +
  facet_wrap(~Model)  + xlab("") 


### Probability of not deviating from +1sd of the mean in protected individuals (test cohort)
dfstore_ind2 <- dfstore_ind[dfstore_ind$Trueclass == "X0",]

dfstore_ind2$notbeX2 <- 1-dfstore_ind2$probX2

dfstore_ind2$OmicsModeling <- "Methylation CpG signature"
dfstore_ind2[dfstore_ind2$Omics== "Transcriptomics" &
               dfstore_ind2$Supervision== "SSAE","OmicsModeling"] <- "Transcriptomics (SSAE)"
dfstore_ind2[dfstore_ind2$Omics== "Transcriptomics" &
               dfstore_ind2$Supervision== "USAE","OmicsModeling"] <- "Transcriptomics (USAE)"
dfstore_ind2[dfstore_ind2$Omics== "Metabolomics" &
               dfstore_ind2$Supervision== "SSAE","OmicsModeling"] <- "Metabolomics (SSAE)"
dfstore_ind2[dfstore_ind2$Omics== "Metabolomics" &
               dfstore_ind2$Supervision== "USAE","OmicsModeling"] <- "Metabolomics (USAE)"
dfstore_ind2[dfstore_ind2$Omics== "Multi-omics" &
               dfstore_ind2$Supervision== "None","OmicsModeling"] <- "Multi-omics"


dfstore_ind2$Model <- as.factor(dfstore_ind2$Model)
ggplot(dfstore_ind2, aes(target, notbeX2, fill= OmicsModeling)) +
  geom_boxplot(aes(group_by=Omics),outlier.shape=" ",alpha=0.9,col="black") +
  theme_classic() +
  scale_fill_manual(values=c("#481568FF" ,"#453781FF","darkgrey","#FDE725FF",
                             "#287D8EFF","#1F968BFF")) +
  ylab("Probability of not deviating from +1sd of the mean in protected individuals (test cohort)") +
  scale_y_continuous(limits = c(0.7, 1)) +
  coord_flip() +
  facet_wrap(~Model)  + xlab("") 



## below are plots used as supplementary figures 

### SBP
tempff <- tempf[tempf$target=="SBP",]
level_order <- factor(tempff$Omics, level = c('Methylation', 'Metabolomics (USAE)', 
                                             'Metabolomics (SSAE)','Transcriptomics (USAE)',
                                             'Transcriptomics (SSAE)',"Multi-omics"))
levels(tempff$`1-sd class F1 (% macro-F1)`) <- c("F1 score specific to unbalanced 1-sd classes (>1sd or <1sd)","F1 score specific to non-deviant 1-sd class (within 1sd)")
tempff$`The macro-F1 score is broken down into (%):` <- tempff$`1-sd class F1 (% macro-F1)`
p <- ggplot(data=tempff, 
       aes(x=level_order, y=F1)) +
  geom_bar(stat="identity", aes(fill=`The macro-F1 score is broken down into (%):`),alpha=0.9)+
  theme_bw() + scale_fill_manual(values=c("firebrick3","grey25"),
                                 guide = guide_legend()) +
  theme(legend.position="bottom",legend.direction="vertical")+
  ylab("Overall performance (macro-F1)") +
  facet_wrap(~Model) + 
  theme(axis.text.x = element_text(angle = 0, hjust=1,size=9),
        axis.text.y = element_text(angle = 0, hjust=1,size=9),
        legend.title=element_text(size=10)) +
  xlab("") + coord_flip() 

p
#ggsave(filename = "plot_SBP.tiff",plot=p, width = 28, height = 14, device='tiff', dpi=900, units = "cm")

### DBP
tempff <- tempf[tempf$target=="DBP",]
level_order <- factor(tempff$Omics, level = c('Methylation', 'Metabolomics (USAE)', 
                                              'Metabolomics (SSAE)','Transcriptomics (USAE)',
                                              'Transcriptomics (SSAE)',"Multi-omics"))
levels(tempff$`1-sd class F1 (% macro-F1)`) <- c("F1 score specific to unbalanced 1-sd classes (>1sd or <1sd)","F1 score specific to non-deviant 1-sd class (within 1sd)")
tempff$`The macro-F1 score is broken down into (%):` <- tempff$`1-sd class F1 (% macro-F1)`
p <- ggplot(data=tempff, 
            aes(x=level_order, y=F1)) +
  geom_bar(stat="identity", aes(fill=`The macro-F1 score is broken down into (%):`),alpha=0.9)+
  theme_bw() + scale_fill_manual(values=c("firebrick3","grey25"),
                                 guide = guide_legend()) +
  theme(legend.position="bottom",legend.direction="vertical")+
  ylab("Overall performance (macro-F1)") +
  facet_wrap(~Model) + 
  theme(axis.text.x = element_text(angle = 0, hjust=1,size=9),
        axis.text.y = element_text(angle = 0, hjust=1,size=9),
        legend.title=element_text(size=10)) +
  xlab("") + coord_flip() 

p
#ggsave(filename = "plot_DBP.tiff",plot=p, width = 28, height = 14, device='tiff', dpi=900, units = "cm")

### EERATIO
tempff <- tempf[tempf$target=="EERATIO",]
level_order <- factor(tempff$Omics, level = c('Methylation', 'Metabolomics (USAE)', 
                                              'Metabolomics (SSAE)','Transcriptomics (USAE)',
                                              'Transcriptomics (SSAE)',"Multi-omics"))
levels(tempff$`1-sd class F1 (% macro-F1)`) <- c("F1 score specific to unbalanced 1-sd classes (>1sd or <1sd)","F1 score specific to non-deviant 1-sd class (within 1sd)")
tempff$`The macro-F1 score is broken down into (%):` <- tempff$`1-sd class F1 (% macro-F1)`
p <- ggplot(data=tempff, 
            aes(x=level_order, y=F1)) +
  geom_bar(stat="identity", aes(fill=`The macro-F1 score is broken down into (%):`),alpha=0.9)+
  theme_bw() + scale_fill_manual(values=c("firebrick3","grey25"),
                                 guide = guide_legend()) +
  theme(legend.position="bottom",legend.direction="vertical")+
  ylab("Overall performance (macro-F1)") +
  facet_wrap(~Model) + 
  theme(axis.text.x = element_text(angle = 0, hjust=1,size=9),
        axis.text.y = element_text(angle = 0, hjust=1,size=9),
        legend.title=element_text(size=10)) +
  xlab("") + coord_flip() 

p
#ggsave(filename = "plot_EERATIO.tiff",plot=p, width = 28, height = 14, device='tiff', dpi=900, units = "cm")

### EARATIO
tempff <- tempf[tempf$target=="EARATIO",]
level_order <- factor(tempff$Omics, level = c('Methylation', 'Metabolomics (USAE)', 
                                              'Metabolomics (SSAE)','Transcriptomics (USAE)',
                                              'Transcriptomics (SSAE)',"Multi-omics"))
levels(tempff$`1-sd class F1 (% macro-F1)`) <- c("F1 score specific to unbalanced 1-sd classes (>1sd or <1sd)","F1 score specific to non-deviant 1-sd class (within 1sd)")
tempff$`The macro-F1 score is broken down into (%):` <- tempff$`1-sd class F1 (% macro-F1)`
p <- ggplot(data=tempff, 
            aes(x=level_order, y=F1)) +
  geom_bar(stat="identity", aes(fill=`The macro-F1 score is broken down into (%):`),alpha=0.9)+
  theme_bw() + scale_fill_manual(values=c("firebrick3","grey25"),
                                 guide = guide_legend()) +
  theme(legend.position="bottom",legend.direction="vertical")+
  ylab("Overall performance (macro-F1)") +
  facet_wrap(~Model) + 
  theme(axis.text.x = element_text(angle = 0, hjust=1,size=9),
        axis.text.y = element_text(angle = 0, hjust=1,size=9),
        legend.title=element_text(size=10)) +
  xlab("") + coord_flip() 

p
#ggsave(filename = "plot_EARATIO.tiff",plot=p, width = 28, height = 14, device='tiff', dpi=900, units = "cm")


### LAVI
tempff <- tempf[tempf$target=="LAVI",]
level_order <- factor(tempff$Omics, level = c('Methylation', 'Metabolomics (USAE)', 
                                              'Metabolomics (SSAE)','Transcriptomics (USAE)',
                                              'Transcriptomics (SSAE)',"Multi-omics"))
levels(tempff$`1-sd class F1 (% macro-F1)`) <- c("F1 score specific to unbalanced 1-sd classes (>1sd or <1sd)","F1 score specific to non-deviant 1-sd class (within 1sd)")
tempff$`The macro-F1 score is broken down into (%):` <- tempff$`1-sd class F1 (% macro-F1)`
p <- ggplot(data=tempff, 
            aes(x=level_order, y=F1)) +
  geom_bar(stat="identity", aes(fill=`The macro-F1 score is broken down into (%):`),alpha=0.9)+
  theme_bw() + scale_fill_manual(values=c("firebrick3","grey25"),
                                 guide = guide_legend()) +
  theme(legend.position="bottom",legend.direction="vertical")+
  ylab("Overall performance (macro-F1)") +
  facet_wrap(~Model) + 
  theme(axis.text.x = element_text(angle = 0, hjust=1,size=9),
        axis.text.y = element_text(angle = 0, hjust=1,size=9),
        legend.title=element_text(size=10)) +
  xlab("") + coord_flip() 

p
#ggsave(filename = "plot_LAVI.tiff",plot=p, width = 28, height = 14, device='tiff', dpi=900, units = "cm")

