## help can be found at: https://tensorflow.rstudio.com/guides/keras/transfer_learning

## load data
setwd("xxx")
Transcripto_FTC <- read.csv("Transcripto_FTC.txt", sep="")
load("xxxY_FIMM.RData")

### load pre-trained SSAE in YFS
library(keras)
library(tensorflow)
model <- load_model_hdf5("transcripto_mod_ss.h5")
summary(model)
pred <- predict(model,as.matrix(Transcripto_FTC))

# test reconstruction
cor(pred[[1]][,5],Transcripto_FTC[,5])

pred.BP <- pred[[2]][,1:2]
cor(pred.BP[,1],Y[,1])
cor(pred.BP[,2],Y[,2])

# take the encoder part of SSAE
encoder <- model$layers[[2]]
pred.encoder <- predict(encoder,as.matrix(Transcripto_FTC))

second = keras_model_sequential()
second <- encoder$output %>% layer_dense(2, input_shape=(6), activation='relu')

summary(second)


################################

## add custom layers
predictions <- encoder$output %>%
  layer_dropout(1/10) %>%
  layer_dense(100, trainable = T) %>%
  layer_activation("relu", trainable = T) %>%
  layer_dense(1, trainable=T) %>%    
  layer_activation("sigmoid", trainable=T)
# this is the model to train
model_f <- keras_model(inputs = encoder$input, outputs = predictions)

# we force pre-trained layers to be untrainable
for (layer in encoder$layers){
  layer$trainable <- F
}

model_f %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.01)
)



X.train <- as.matrix(scale(Transcripto_FTC))
### create hypertension categories
X.train[X.train=="NaN"] <- 0
Y.train <- Y[,2]
Y.sys <- Y.train
Y.sys[Y.train>=160] <- 1
Y.sys[Y.train<160] <- 0
table(Y.sys)
Y.sys <- as.numeric(Y.sys)

# fit model
history <- model_f %>% fit(as.matrix(X.train),Y.sys, epochs=100)

pred<- predict(model_f,X.test)
plot(pred,Ysys.test)
# get ROC
library(pROC)
test_roc = roc(Ysys.test ~ as.numeric(pred), plot = TRUE, print.auc = TRUE)


######
# this is a clone model
model_compare <- clone_model(model_f)


model_compare %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.001)  ## play with the learning rate

)
history <- model_compare %>% fit(as.matrix(X.train),Y.sys, epochs=100)
pred<- predict(model_compare,X.train)
plot(pred,Y.sys)
library(pROC)
test_roc = roc(Y.sys ~ as.numeric(pred), plot = TRUE, print.auc = TRUE)



#### what about using a smaller training samples?

i <- sample(1:310,replace = F)[1:100]
X.train <- Transcripto_FTC[i,]
Ysys.train <- Y.sys[i]
X.test <- Transcripto_FTC[-i,]
Ysys.test <- Y.sys[-i]

## scaling
for(k in 1:ncol(X.train)){
  m <- mean(X.train[,k]) ; sd <- sd(X.train[,k])
  X.train[,k] <- (X.train[,k]-m)/sd
  X.test[,k] <- (X.test[,k]-m)/sd
}
X.train[X.train=="NaN"] <- 0 #if any
X.test[X.test=="NaN"] <- 0

X.train <- as.matrix(X.train)
X.test <- as.matrix(X.test)

##model :
predictions <- encoder$output %>%
  layer_dropout(1/6) %>%
  layer_dense(6, trainable = T) %>%
  layer_activation("relu", trainable = T) %>%
  layer_dense(1, trainable=T) %>%    
  layer_activation("sigmoid", trainable=T)
model_f <- keras_model(inputs = encoder$input, outputs = predictions)
for (layer in encoder$layers){
  layer$trainable <- T
}
model_f %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.01),  
  metrics = 'accuracy'
)
history <- model_f %>% fit(as.matrix(X.train),Ysys.train, epochs=100)

pred<- predict(model_f,X.test)
plot(pred,Ysys.test)
library(pROC)
test_roc = roc(Ysys.test ~ as.numeric(pred), plot = TRUE, print.auc = TRUE)


## 
predictions <- encoder$output %>%
  layer_dropout(1/6) %>%
  layer_dense(6, trainable = T) %>%
  layer_activation("relu", trainable = T) %>%
  layer_dense(1, trainable=T) %>%    
  layer_activation("sigmoid", trainable=T)
model_f <- keras_model(inputs = encoder$input, outputs = predictions)
for (layer in encoder$layers){
  layer$trainable <- T
}
model_f %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.01),  ## play with the learning rate
  metrics = 'accuracy'
)
model_compare <- clone_model(model_f)


model_compare %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.001), 
  metrics = 'accuracy'
)
history <- model_compare %>% fit(as.matrix(X.train),Ysys.train, epochs=100)
pred<- predict(model_compare,X.test)
plot(pred,Ysys.test)
library(pROC)
test_roc = roc(Ysys.test ~ as.numeric(pred), plot = TRUE, print.auc = TRUE)


### what about adding age + BMI + sex ????
# we load the data
load("X_FIMM.RData")
C <- X$Clinical_PRS[,c(2,22,24)]
i <- sample(1:310,replace = F)[1:200]
X.train <- Transcripto_FTC[i,]
Ysys.train <- Y.sys[i]
X.test <- Transcripto_FTC[-i,]
Ysys.test <- Y.sys[-i]
C.train <- C[i,]
C.test <- C[-i,]


## scaling
for(k in 1:ncol(X.train)){
  m <- mean(X.train[,k]) ; sd <- sd(X.train[,k])
  X.train[,k] <- (X.train[,k]-m)/sd
  X.test[,k] <- (X.test[,k]-m)/sd
}
X.train[X.train=="NaN"] <- 0
X.test[X.test=="NaN"] <- 0

X.train <- as.matrix(X.train)
X.test <- as.matrix(X.test)

for(k in 2:ncol(C.train)){
  m <- mean(C.train[,k]) ; sd <- sd(C.train[,k])
  C.train[,k] <- (C.train[,k]-m)/sd
  C.test[,k] <- (C.test[,k]-m)/sd
}

inpclin <- layer_input(shape=c(3))
newinc <- layer_concatenate(list(encoder$output,inpclin))

predictions <- encoder$output %>%
  layer_dropout(1/10) %>%
  layer_dense(100, trainable = T) %>%
  layer_activation("relu", trainable = T) 

predictions2 <- layer_concatenate(list(predictions,inpclin)) %>%
  layer_dense(1, trainable=T) %>%    
  layer_activation("sigmoid", trainable=T)

model_f <- keras_model(inputs = list(encoder$input,inpclin), outputs = predictions2)

for (layer in encoder$layers){
  layer$trainable <- T
}
model_f %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.01),  ## play with the learning rate
  metrics = 'accuracy'
)

history <- model_f %>% fit(list(as.matrix(X.train), as.matrix(C.train)),
                           Ysys.train, epochs=100)

pred<- predict(model_f,list(as.matrix(X.test),as.matrix(C.test)))
plot(pred,Ysys.test)
library(pROC)
test_roc = roc(Ysys.test ~ as.numeric(pred), plot = TRUE, print.auc = TRUE)

##############################################################
############### get T-sne representation for a layer

### Test sample
layer_name <- 'Mylayer'
intermediate_layer_model <- keras_model(inputs = model_f$input,
                                        outputs = get_layer(model_f, layer_name)$output)
intermediate_output <- predict(intermediate_layer_model, list(as.matrix(X.test),as.matrix(C.test)))

require(Rtsne)

mytsne <- Rtsne(intermediate_output,perplexity = 10)
dftsne <- as.data.frame(mytsne$Y)
dftsne$Status <- as.factor(Ysys.test)


require(ggplot2)
ggplot(dftsne,aes(V1,V2,color=Status)) + geom_point(size=4,shape=19) +
  theme_void() + scale_color_manual(values=c("black","red")) +
  ggtitle("") 

# Train sample
layer_name <- 'Mylayer'
intermediate_layer_model <- keras_model(inputs = model_f$input,
                                        outputs = get_layer(model_f, layer_name)$output)
intermediate_output <- predict(intermediate_layer_model, list(as.matrix(X.train),as.matrix(C.train)))

require(Rtsne)

set.seed(1871)
mytsne <- Rtsne(intermediate_output,perplexity = 20)
dftsne <- as.data.frame(mytsne$Y)
dftsne$Status <- as.factor(Ysys.train)
levels(dftsne$Status) <- c("Control","Hypertensive")

# 

require(ggplot2)
dftsne <- dftsne[order(dftsne$Status),]
ggplot(dftsne,aes(V1,V2,color=Status)) + geom_point(size=4,shape=19,alpha=0.7) +
  theme_test() + scale_color_manual(values=c("lightgrey","darkred")) +
  ggtitle("2D t-SNE representation | Transcriptomics") + 
  theme(legend.position="right",
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=13)) +
  xlab("dimension 1") + ylab("dimension 2") +
  annotate("text", x=-12, y=-9, label= "Training proportion: 80%",size=3.5) +
  annotate("text", x=-8.4, y=-11.5, label= "Last MLP-layer before concatenation", size=3.5) +
  xlim(-18,18) 

ggsave(filename = "plottsne.tiff", device='tiff', dpi=600)

