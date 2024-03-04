# we build a function to automatize the use of ML classifiers
# ss sufixes: SSAE. us suffixes: USAE.

# xtarget: CVD outcome
# xmethod: classifier to be used
# xdata.ss/us: training set
# xdatatest.ss/us: testing set
# xtruetest: "true" classes in test sample (N=250)
# xdfstore_ind/xdfstore_mod: datasets to fill in (individual predictions or model performance)

trainfunction <- function(xtarget,xmethod,xdata.ss,xdata.us,xdatatest.ss,
                          xdatatest.us,xtruetest,xdfstore_ind,xdfstore_mod){
  ### SSAE first
  mod.ss <- caret::train(target~., data=xdata.ss,
                         method=xmethod,
                         metric = "F1",
                         trControl = train.control)
  
  pred.ss <- predict(mod.ss, newdata = xdatatest.ss,type = "raw")
  xdfstore_mod <- rbind(xdfstore_mod,
                        c(xtarget,xmethod,"SSAE",
                          f1class(table(mod.ss$pred$pred,mod.ss$pred$obs)),
                          f1class(table(pred.ss,xtruetest))))
  
  pred.ss <- predict(mod.ss, newdata = xdatatest.ss,type = "prob")
  
  xdfstore_ind <- rbind(xdfstore_ind,
                        data.frame(
                          target=rep(xtarget,249),
                          Model=rep(xmethod,249),
                          Supervision=rep("SSAE",249),
                          Trueclass=xtruetest,
                          probX0=pred.ss[,1],
                          probX1=pred.ss[,2],
                          probX2=pred.ss[,3]))
  
  pred.ss <- predict(mod.ss, newdata = xdata.ss,type = "prob")
  
  save_train.ss <- data.frame( 
    target=rep(xtarget,1000),
    Model=rep(xmethod,1000),
    Supervision=rep("SSAE",1000),
    Trueclass=xdata.ss$target,
    probX0=pred.ss[,1],
    probX1=pred.ss[,2],
    probX2=pred.ss[,3])
  
  ### USAE then
  
  mod.us <- caret::train(target~., data=xdata.us,
                         method=xmethod,
                         metric = "F1",
                         trControl = train.control)
  
  
  pred.us <- predict(mod.us, newdata = xdatatest.us,type = "raw")
  xdfstore_mod <- rbind(xdfstore_mod,
                        c(xtarget,xmethod,"USAE",
                          f1class(table(mod.us$pred$pred,mod.us$pred$obs)),
                          f1class(table(pred.us,xtruetest))))
  
  pred.us <- predict(mod.us, newdata = xdatatest.us,type = "prob") 
  #type prob is optional (one could use classes, but for VCDN we need the probabilities of belonging in a class)
  
  xdfstore_ind <- rbind(xdfstore_ind,
                        data.frame(
                          target=rep(xtarget,249),
                          Model=rep(xmethod,249),
                          Supervision=rep("USAE",249),
                          Trueclass=xtruetest,
                          probX0=pred.us[,1],
                          probX1=pred.us[,2],
                          probX2=pred.us[,3]))
  
  pred.us <- predict(mod.us, newdata = xdata.us,type = "prob")
  
  save_train.us <- data.frame( 
    target=rep(xtarget,1000),
    Model=rep(xmethod,1000),
    Supervision=rep("USAE",1000),
    Trueclass=xdata.us$target,
    probX0=pred.us[,1],
    probX1=pred.us[,2],
    probX2=pred.us[,3])
  
  return(list(xdfstore_ind,xdfstore_mod,save_train.ss,save_train.us))
}

#########################
## preparing a dataset to store the results for each individual 

dfstore_ind <- data.frame(matrix(NA,ncol=7))
colnames(dfstore_ind) <- c("target","Model","Supervision","Trueclass","probX0","probX1","probX2")

dfstore_mod <- data.frame(matrix(NA,ncol=9))
colnames(dfstore_mod) <- c("target","Model","Supervision",
                           "F1trainX0","F1trainX1","F1trainX2",
                           "F1testX0","F1testX1","F1testX2")

set.seed(1871)

#######################################################
#### We take the example of systolic blood pressure (SBP)
### one can change the outcome and re-run models.

target <- Y.train[,6]
truetest <- Y.test[,6]

XT.trss <- as.data.frame(XT.trss)
XT.trus <- as.data.frame(XT.trus)
XMb.trss <- as.data.frame(XMb.trss)
XMb.trus <- as.data.frame(XMb.trus)
XMe.tr <- as.data.frame(XMe.tr)

XT.trss$target <- target
XT.trus$target <- target
XMb.trss$target <- target
XMb.trus$target <- target
XMe.tr$target <- target

########################

#rf
temp <- trainfunction("SBP","rf",
                      xdata.ss=XT.trss,
                      xdata.us = XT.trus,
                      xdatatest.ss = XT.tess,
                      xdatatest.us = XT.teus,
                      xtruetest = truetest,
                      xdfstore_ind = dfstore_ind,
                      xdfstore_mod = dfstore_mod)
dfstore_ind <- as.data.frame(temp[[1]])
dfstore_mod <- as.data.frame(temp[[2]])
dfstore_trainingind <- rbind(as.data.frame(temp[[3]]),
                             as.data.frame(temp[[4]]))

# lda
temp <- trainfunction("SBP","lda",
                      xdata.ss=XT.trss,
                      xdata.us = XT.trus,
                      xdatatest.ss = XT.tess,
                      xdatatest.us = XT.teus,
                      xtruetest = truetest,
                      xdfstore_ind = dfstore_ind,
                      xdfstore_mod = dfstore_mod)
dfstore_ind <- as.data.frame(temp[[1]])
dfstore_mod <- as.data.frame(temp[[2]])
dfstore_trainingind <- rbind(dfstore_trainingind,
                             as.data.frame(temp[[3]]),
                             as.data.frame(temp[[4]]))



# gbm
temp <- trainfunction("SBP","gbm",
                      xdata.ss=XT.trss,
                      xdata.us = XT.trus,
                      xdatatest.ss = XT.tess,
                      xdatatest.us = XT.teus,
                      xtruetest = truetest,
                      xdfstore_ind = dfstore_ind,
                      xdfstore_mod = dfstore_mod)
dfstore_ind <- as.data.frame(temp[[1]])
dfstore_mod <- as.data.frame(temp[[2]])
dfstore_trainingind <- rbind(dfstore_trainingind,
                             as.data.frame(temp[[3]]),
                             as.data.frame(temp[[4]]))

# svm linear 
temp <- trainfunction("SBP","svmLinear",
                      xdata.ss=XT.trss,
                      xdata.us = XT.trus,
                      xdatatest.ss = XT.tess,
                      xdatatest.us = XT.teus,
                      xtruetest = truetest,
                      xdfstore_ind = dfstore_ind,
                      xdfstore_mod = dfstore_mod)
dfstore_ind <- as.data.frame(temp[[1]])
dfstore_mod <- as.data.frame(temp[[2]])
dfstore_trainingind <- rbind(dfstore_trainingind,
                             as.data.frame(temp[[3]]),
                             as.data.frame(temp[[4]]))

# svm radial 
temp <- trainfunction("SBP","svmRadial",
                      xdata.ss=XT.trss,
                      xdata.us = XT.trus,
                      xdatatest.ss = XT.tess,
                      xdatatest.us = XT.teus,
                      xtruetest = truetest,
                      xdfstore_ind = dfstore_ind,
                      xdfstore_mod = dfstore_mod)
dfstore_ind <- as.data.frame(temp[[1]])
dfstore_mod <- as.data.frame(temp[[2]])
dfstore_trainingind <- rbind(dfstore_trainingind,
                             as.data.frame(temp[[3]]),
                             as.data.frame(temp[[4]]))

# perceptron 
temp <- trainfunction("SBP","mlp",
                      xdata.ss=XT.trss,
                      xdata.us = XT.trus,
                      xdatatest.ss = XT.tess,
                      xdatatest.us = XT.teus,
                      xtruetest = truetest,
                      xdfstore_ind = dfstore_ind,
                      xdfstore_mod = dfstore_mod)
dfstore_ind <- as.data.frame(temp[[1]])
dfstore_mod <- as.data.frame(temp[[2]])
dfstore_trainingind <- rbind(dfstore_trainingind,
                             as.data.frame(temp[[3]]),
                             as.data.frame(temp[[4]]))

