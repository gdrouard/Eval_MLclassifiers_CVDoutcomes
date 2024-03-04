dfstore_mod <- read.csv("xxx/dfstore_mod.txt", sep="")

## remove F1-score from training set
colnames(dfstore_mod)
dfstore_mod <- dfstore_mod[,-c(4,5,6)]

## add multiomic predictions
summary_Multiomics_pred <- read.csv("xxx/summary_Multiomics_pred.txt", sep="")
colnames(summary_Multiomics_pred)
summary_Multiomics_pred <- summary_Multiomics_pred[,c(1,2,4,5,6)]
summary_Multiomics_pred$Supervison <- "SSAE"
summary_Multiomics_pred <- summary_Multiomics_pred[,c(1,2,6,3,4,5)]

# rbind datasets
colnames(dfstore_mod)
colnames(summary_Multiomics_pred)
colnames(summary_Multiomics_pred) <- colnames(dfstore_mod)
dfstore_mod <- as.data.frame(rbind(dfstore_mod,summary_Multiomics_pred))


dfstore_mod[dfstore_mod=="NaN"] <- 0
dfstore_mod$Model <- as.factor(dfstore_mod$Model)
levels(dfstore_mod$Model) <- c("Gradient Boosting Machine","Linear Discriminant Analysis",
                                "Multi-layer Perceptron","Random forest","Linear SVM","Radial SVM")
dfstore_mod$Supervision <- as.factor(dfstore_mod$Supervision)
levels(dfstore_mod$Supervision) <- c("None (feature selection - DNA methylation)","Semi-supervised AE",
                               "Unsupervised AE")

dfstore_mod$Variable <- paste0(dfstore_mod$Model,"|",dfstore_mod$Supervision)

dfstore_mod$target <- as.factor(dfstore_mod$target)
levels(dfstore_mod$target) <- c("Diastolic BP","E/A ratio","E/e' ratio","LAVI","Systolic BP")

### splitting dataset before long format (for each 1-sd class)
colnames(dfstore_mod)
dfX0 <- dfstore_mod[,c(1,2,3,7,8,4)]
dfX1 <- dfstore_mod[,c(1,2,3,7,8,5)]
dfX2 <- dfstore_mod[,c(1,2,3,7,8,6)]

colnames(dfX0) <- c("target","Model","Supervision","Omic","Variable","F1")
colnames(dfX1) <- c("target","Model","Supervision","Omic","Variable","F1")
colnames(dfX2) <- c("target","Model","Supervision","Omic","Variable","F1")

dfX0$`1-sd classes` <- "<1sd from mean"
dfX1$`1-sd classes` <- "within 1sd from mean"
dfX2$`1-sd classes` <- ">1sd from mean"

#merging
df <- rbind(dfX0,dfX1,dfX2)

head(df)

###################################
###########   PLOTS  ##############
###################################


setwd("xxx") #where to save figures

require(ggplot2)

ggplot(df[df$Supervision=="Semi-supervised AE" & df$Omic=="Transcriptomics",]) +
  geom_bar( aes(x=Model, y=F1, fill= `1-sd classes`), stat="identity", alpha=0.7,position = position_dodge()) +
  facet_grid(~target) + coord_flip() + theme_bw()

ggplot(df[(df$Supervision=="Semi-supervised AE" | 
             df$Supervision=="None (feature selection - DNA methylation)") & 
            df$`1-sd classes`=="<1sd from mean",]) +
  geom_bar( aes(x=target, y=F1, fill= Omic), stat="identity", alpha=0.7,position = position_dodge()) +
  facet_wrap(~Model) + coord_flip() + theme_bw() + xlab(" ") +
  ggtitle("Evaluation of predictions of being <1sd from the mean")+
  ylab("class-specific F1 score")+
  scale_fill_manual(values=c("darkcyan","darkorchid4","grey10","firebrick3"))+
  theme(legend.position = "bottom",legend.title=element_blank(),
        legend.text = element_text(size=12, face="bold"))

ggsave(filename = "rev2_fig1.png" ,width = 22, height = 13, device='png', dpi=500,units = "cm")

ggplot(df[(df$Supervision=="Semi-supervised AE" | 
             df$Supervision=="None (feature selection - DNA methylation)") & 
            df$`1-sd classes`==">1sd from mean",]) +
  geom_bar( aes(x=target, y=F1, fill= Omic), stat="identity", alpha=0.7,position = position_dodge()) +
  facet_wrap(~Model) + coord_flip() + theme_bw() + xlab(" ") +
  ggtitle("Evaluation of predictions of being >1sd from the mean") +
  ylab("class-specific F1 score") +
  scale_fill_manual(values=c("darkcyan","darkorchid4","grey10","firebrick3"))+
  theme(legend.position = "bottom",legend.title=element_blank(),
        legend.text = element_text(size=12, face="bold"))

ggsave(filename = "rev2_fig2.png" ,width = 22, height = 13, device='png', dpi=500,units = "cm")

