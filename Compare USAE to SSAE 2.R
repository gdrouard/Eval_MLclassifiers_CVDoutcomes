
# LOAD SUMMARY STATS
#summary stats single-omics:
dfstore_mod <- read.csv("xxx/dfstore_mod.txt", sep="") 


#####################################################
# Investigate whether SSAE outperforms USAE   #######
#####################################################

head(dfstore_mod)

### In both transcriptomic and metabolomic levels:
temp <- dfstore_mod[dfstore_mod$Omics %in% c("Transcriptomics","Metabolomics"),]
head(temp)
tempssae <- temp[temp$Supervision=="SSAE",]
tempusae <- temp[temp$Supervision=="USAE",]
mean(tempssae$F1);mean(tempusae$F1)
hist(tempssae$F1-tempusae$F1)
quantile(tempssae$F1-tempusae$F1,probs=myprobs)


### At single-omic level:
temp <- dfstore_mod[dfstore_mod$Omics == "Metabolomics",] #change by Transcriptomics if needed
head(temp)
tempssae <- temp[temp$Supervision=="SSAE",]
tempusae <- temp[temp$Supervision=="USAE",]
mean(tempssae$F1);mean(tempusae$F1)
myprobs<-seq(0,1,length=nrow(tempssae))
hist(tempssae$F1-tempusae$F1)
quantile(tempssae$F1-tempusae$F1,probs=myprobs)
mean(tempssae$F1-tempusae$F1)
quantile(tempssae$F1testX2-tempusae$F1testX2,probs=myprobs)
mean(tempssae$F1testX2-tempusae$F1testX2)
quantile(tempssae$F1testX0-tempusae$F1testX0,probs=myprobs)
mean(tempssae$F1testX0-tempusae$F1testX0)

### Proportion of SSAE > USAE (%)
length(which(tempssae$F1-tempusae$F1==0))/nrow(tempssae)*100 # equal
length(which(tempssae$F1-tempusae$F1>0))/nrow(tempssae)*100 #strictly superior
length(which(tempssae$F1testX0-tempusae$F1testX0>0))/nrow(tempssae)*100
length(which(tempssae$F1testX2-tempusae$F1testX2>0))/nrow(tempssae)*100


