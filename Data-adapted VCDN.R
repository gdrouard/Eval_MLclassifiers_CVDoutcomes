### training for VCDN

library(keras)
library(tensorflow)

# function. Target is CVD outcome, MLmodel is which ML classifier to be used.
VCDN_training_fun <- function(target,MLmodel){
  
  ### know how to take target information:
  if(target=="SBP"){varindex <- c(1:3)}
  if(target=="DBP"){varindex <- c(4:6)}
  if(target=="EERATIO"){varindex <- c(7:9)}
  if(target=="EARATIO"){varindex <- c(10:12)}
  if(target=="LAVI"){varindex <- c(13:15)}
  
  ### Prepare data
  
  # dfstore_trainingind: individual-level data
  dftemp <- dfstore_trainingind[dfstore_trainingind$Model == MLmodel &
                                  dfstore_trainingind$Supervision %in% c("SSAE","None"),]
  temp <- dftemp[dftemp$Omics == "Transcriptomics",]
  # temp is in long format, with order target variables. We want it to be wide format
  # which can be done manually or with e.g. dplyr. Below is manual version.
  BPpred.T <- cbind(
    temp[1:1000,5:7],
    temp[1001:2000,5:7],
    temp[2001:3000,5:7],
    temp[3001:4000,5:7],
    temp[4001:5000,5:7]
  )
  colnames(BPpred.T) <- c("X0.sbp","X1.sbp","X2.sbp",
                          "X0.dbp","X1.dbp","X2.dbp",
                          "X0.eeratio","X1.eeratio","X2.eeratio",
                          "X0.earatio","X1.earatio","X2.earatio",
                          "X0.lavi","X1.lavi","X2.lavi")
  
  # similar process for metabolomic and methylation data
  temp <- dftemp[dftemp$Omics == "Metabolomics",]
  BPpred.Mb <- cbind(
    temp[1:1000,5:7],
    temp[1001:2000,5:7],
    temp[2001:3000,5:7],
    temp[3001:4000,5:7],
    temp[4001:5000,5:7]
  )
  colnames(BPpred.Mb) <- colnames(BPpred.T)
  BPpred.Me <- cbind(
    temp[1:1000,5:7],
    temp[1001:2000,5:7],
    temp[2001:3000,5:7],
    temp[3001:4000,5:7],
    temp[4001:5000,5:7]
  )
  colnames(BPpred.Me) <- colnames(BPpred.T)
  
  # prepare true class distribution of individuals
  # we take methylation data as arbitrary data asindividuals are ordered the same way for all omics
  temp <- dftemp[dftemp$Omics == "Methylation",]
  Truedistribution <- cbind(
    temp[1:1000,"Trueclass"],
    temp[1001:2000,"Trueclass"],
    temp[2001:3000,"Trueclass"],
    temp[3001:4000,"Trueclass"],
    temp[4001:5000,"Trueclass"]
  )
  
  ### encode Y matrix into one-hot encoded label matrix
  TruedistribOH <- matrix(0,nrow=1000,ncol=15)
  for(k in 1:1000){
    colk=1
    for(k2 in 1:5){
      if(Truedistribution[k,k2]=="X0"){
        TruedistribOH[k,colk] <- 1
      }
      if(Truedistribution[k,k2]=="X1"){
        TruedistribOH[k,colk+1] <- 1
      }
      if(Truedistribution[k,k2]=="X2"){
        TruedistribOH[k,colk+2] <- 1
      }
      colk=colk+3
    }
  }
  
  ### Prepare resized tensor. Later on, we stored everything into a vector.
  CrossMat <- matrix(NA,nrow=1000,ncol=3*3*3)
  BPpred.T <- BPpred.T[,varindex] #take target info
  BPpred.Mb <- BPpred.Mb[,varindex] #take target info
  BPpred.Me <- BPpred.Me[,varindex] #take target info
  for(nb in 1:1000){
    tensor_resized <- rep(NA,3*3*3)
    for(i in 1:3){
      for(j in 1:3){
        for(k in 1:3){
          index <- (i-1)*9 + (j-1)*3 + k
          ci <- BPpred.T[nb,i] * BPpred.Mb[nb,j] * BPpred.Me[nb,k]
          tensor_resized[index] <- ci
        }
      }
    }
    CrossMat[nb,] <- tensor_resized
  }
  
  ### scale Crossmat and keep the mean/sd info for testing
  msd <- as.data.frame(matrix(NA,ncol=2,nrow=ncol(CrossMat)))
  for(i in 1:ncol(CrossMat)){
    m <- mean(CrossMat[,i]); sd <- sd(CrossMat[,i])
    CrossMat[,i] <- (CrossMat[,i]-m)/sd
    msd[i,] <- c(m,sd)
  }

  # finally take target info for label distribution (for training)
  TruedistribOH <- TruedistribOH[,varindex]
  
  # VCDN Training
  CrossMat <- as.matrix(CrossMat)
  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = ncol(CrossMat), activation = 'relu', 
                input_shape = c(ncol(CrossMat))) %>%
    layer_dropout(rate=0.5) %>% 
    layer_dense(units=9,activation="relu") %>%
    layer_dense(units = 3, activation = 'softmax')
  summary(model)
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(learning_rate =0.001)
  )
  history <- model %>% fit(
    as.matrix(CrossMat),TruedistribOH, 
    epochs = 30, batch_size = 32, 
    validation_split = 0.2
  )
  
  # predictions from training set
  pred <- predict(object = model, x=as.matrix(CrossMat))
  
  #######
  ### Let's test with test sample !!
  dftemp <- dfstore_ind[dfstore_ind$Model == MLmodel &
                          dfstore_ind$Supervision %in% c("SSAE","None"),]
  temp <- dftemp[dftemp$Omics == "Transcriptomics",]
  BPpred.T <- cbind(
    temp[1:249,5:7],
    temp[250:498,5:7],
    temp[499:747,5:7],
    temp[748:996,5:7],
    temp[997:1245,5:7]
  )
  colnames(BPpred.T) <- colnames(BPpred.T)
  temp <- dftemp[dftemp$Omics == "Metabolomics",]
  BPpred.Mb <- cbind(
    temp[1:249,5:7],
    temp[250:498,5:7],
    temp[499:747,5:7],
    temp[748:996,5:7],
    temp[997:1245,5:7]
  )
  colnames(BPpred.Mb) <- colnames(BPpred.T)
  temp <- dftemp[dftemp$Omics == "epigenetics",]
  BPpred.Me <- cbind(
    temp[1:249,5:7],
    temp[250:498,5:7],
    temp[499:747,5:7],
    temp[748:996,5:7],
    temp[997:1245,5:7]
  )
  colnames(BPpred.Me) <- colnames(BPpred.T)
  
  # prepare true distribution of individuals
  temp <- dftemp[dftemp$Omics == "epigenetics",]
  Truedistribution.test <- cbind(
    temp[1:249,"Trueclass"],
    temp[250:498,"Trueclass"],
    temp[499:747,"Trueclass"],
    temp[748:996,"Trueclass"],
    temp[997:1245,"Trueclass"]
  )
  
  TruedistribOH.test <- matrix(0,nrow=249,ncol=15)
  for(k in 1:249){
    colk=1
    for(k2 in 1:5){
      if(Truedistribution.test[k,k2]=="X0"){
        TruedistribOH.test[k,colk] <- 1
      }
      if(Truedistribution.test[k,k2]=="X1"){
        TruedistribOH.test[k,colk+1] <- 1
      }
      if(Truedistribution.test[k,k2]=="X2"){
        TruedistribOH.test[k,colk+2] <- 1
      }
      colk=colk+3
    }
  }
  
  CrossMat.test <- matrix(NA,nrow=249,ncol=3*3*3)
  BPpred.T <- BPpred.T[,varindex]
  BPpred.Mb <- BPpred.Mb[,varindex]
  BPpred.Me <- BPpred.Me[,varindex]
  
  for(nb in 1:249){
    tensor_resized <- rep(NA,3*3*3)
    for(i in 1:3){
      for(j in 1:3){
        for(k in 1:3){
          index <- (i-1)*9 + (j-1)*3 + k
          ci <- BPpred.T[nb,i] * BPpred.Mb[nb,j] * BPpred.Me[nb,k]
          tensor_resized[index] <- ci
        }
      }
    }
    CrossMat.test[nb,] <- tensor_resized
  }
  
  ### scale Crossmat with the mean/sd info for testing
  for(i in 1:ncol(CrossMat.test)){
    CrossMat.test[,i] <- (CrossMat.test[,i]-msd[i,1])/msd[i,2]
  }
  # keep new predictions
  pred <- predict(object = model, x=as.matrix(CrossMat.test))
  TruedistribOH.test <- TruedistribOH.test[,varindex]
  
  # recode (highest probabilities suggest the class)
  pred2 <- pred
  for(i in 1:249){
    for(k in 1:3){ 
      if(pred[i,k]==max(pred[i,])){pred2[i,k] <- "P1"} 
      if(pred[i,k] != max(pred[i,])){pred2[i,k] <- "P0"}
    }
  }
  pred3 <- pred2[,1]
  for(i in 1:nrow(pred2)){
    if(pred2[i,1]=="P1"){
      pred3[i] <- "X0"
    }
    if(pred2[i,2]=="P1"){
      pred3[i] <- "X1"
    }
    if(pred2[i,3]=="P1"){
      pred3[i] <- "X2"
    }
  }
  
  # true 1-sd class to compare with predictions
  if(target == "SBP"){TD2 <- Truedistribution.test[,1]}
  if(target == "DBP"){TD2 <- Truedistribution.test[,2]}
  if(target == "EERATIO"){TD2 <- Truedistribution.test[,3]}
  if(target == "EARATIO"){TD2 <- Truedistribution.test[,4]}
  if(target == "LAVI"){TD2 <- Truedistribution.test[,5]}
  
  
  # final step: saving results
  t <- table(pred3,TD2)
  if(dim(t)[1]==2 & dim(t)[2]==3){t <- rbind(c(0,0,0),t)}
  if(dim(t)[1]==1 & dim(t)[2]==3){t <- rbind(c(0,0,0),t)}
  f1c <- f1class(t) #f1class function calculates F1 score
  macrof1 <- mean(f1class(t)) 
  
  pred.r <- cbind(rep(target,249),rep(MLmodel,249),
                  rep("Multi-omics",249),TD2,pred)
  colnames(pred.r) <- c("target","Model","Mode","Trueclass","PredX0",
                        "predX1","predX2")
  pred.r <- as.data.frame(pred.r)
  results_ <- list(c(macrof1,f1c),pred.r)
  
  return(results_)
  
}

