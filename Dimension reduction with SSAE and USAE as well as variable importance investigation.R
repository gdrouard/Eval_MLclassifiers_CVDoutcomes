
### building an SSAE

latent_size = 50

enc_input = layer_input(shape = ncol(X.train))
enc_output = enc_input %>% 
  layer_gaussian_noise(0.01) %>%
  layer_dropout(rate=0.5) %>%
  layer_dense(units=150, kernel_constraint = constraint_maxnorm(3)) %>% 
  layer_activation_leaky_relu(alpha=0.1) %>%
  layer_dense(units=latent_size,kernel_constraint = constraint_maxnorm(3)) %>%
  layer_activation_leaky_relu(alpha=0.1)

encoder = keras_model(enc_input, enc_output,name="encoder")

dec_input = layer_input(shape = latent_size)
dec_output = dec_input %>% 
  layer_dense(units=150) %>% 
  layer_activation_leaky_relu(alpha=0.1) %>% 
  layer_dense(units = ncol(X.train), name = "reconstruction_output") %>% 
  layer_activation_leaky_relu(alpha=0.1) 

decoder = keras_model(dec_input, dec_output,name="decoder")
ann_input = layer_input(shape = latent_size)
auxiliary_output <- ann_input %>% 
  layer_dense(units=latent_size) %>%
  layer_activation_leaky_relu(0.1) %>%
  layer_dense(units=ncol(Y.train),activation="linear")
ann <- keras_model(ann_input, auxiliary_output,name="ANN")

aen_input = layer_input(shape = ncol(X.train))
aen_output1 = aen_input %>% 
  encoder() %>% 
  decoder()

aen_output2 = aen_input %>% 
  encoder() %>% 
  ann()

aen = keras_model(aen_input, c(aen_output1,aen_output2),name="Modelfinal")
lr_schedule = optimizer_adam(
  learning_rate=0.001,
  beta_1 =0.9)
aen %>% compile(
  optimizer = lr_schedule,
  loss = 'mse',
  loss_weights = c(0.9,0.1),
  metrics = c("mse")
)
history <- aen %>% fit(X.train,list(X.train,Y.train), epochs=100, batch_size=64, validation_split=0.2) 

#train data 
latentspace.train.ss = encoder %>% predict(X.train)
decoded_space = decoder %>% predict(latentspace.train.ss)
decoded_space[1:5,1:4]
plot(decoded_space[,1],X.train[,1]);abline(a=0,b=1)
plot(decoded_space[,10],X.train[,10]);abline(a=0,b=1)

BPpred.train.ss <- ann %>% predict(latentspace.train.ss)
plot(BPpred.train.ss[,1],Y.train[,1],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
plot(BPpred.train.ss[,2],Y.train[,2],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
cor(BPpred.train.ss[,1],Y.train[,1])
cor(BPpred.train.ss[,2],Y.train[,2])

# test data
latentspace.test.ss = encoder %>% predict(X.test)
decoded_space = decoder %>% predict(latentspace.test.ss)
latentspace.test.ss[1:5,1:4]
plot(decoded_space[,1],X.test[,1]);abline(a=0,b=1)
plot(decoded_space[,10],X.test[,10]);abline(a=0,b=1)

BPpred.test.ss <- ann %>% predict(latentspace.test.ss)
plot(BPpred.test.ss[,1],Y.test[,1],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
plot(BPpred.test.ss[,2],Y.test[,2],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
cor(BPpred.test.ss[,1],Y.test[,1])
cor(BPpred.test.ss[,2],Y.test[,2])
cor(BPpred.test.ss[,3],Y.test[,3])
cor(BPpred.test.ss[,4],Y.test[,4])
cor(BPpred.test.ss[,5],Y.test[,5])

############################
### CROSS-variable IMPORTANCE

W1 <- keras::get_weights(encoder)
W2 <- keras::get_weights(ann)

Wb1 <- W1[[1]] %*% W1[[3]] 
Wb2 <- W2[[1]] %*% W2[[3]] 
#Wb1 <- abs(Wb1); Wb2 <- abs(Wb2) #no
dim(Wb1);dim(Wb2)
Wbcross <- Wb1 %*% Wb2

rownames(Wbcross) <- colnames(X.train)
colnames(Wbcross) <- colnames(Y.train)

### plots
head(Wbcross)
head(Wbcross[order(Wbcross[,1],decreasing=T),1],20)
head(Wbcross[order(Wbcross[,2],decreasing=T),2],20)
head(Wbcross[order(Wbcross[,3],decreasing=T),3],20)
head(Wbcross[order(Wbcross[,4],decreasing=T),4],20)
head(Wbcross[order(Wbcross[,5],decreasing=T),5],20)

# genes to look at
keyrepgenes <- c("MYADM","MYADM.1","CD97","TPPP3","TIPARP","SLC31A2",
                 "LMNA","F12","CRIP1","TSC22D3","CEBPA",
                 "SLC31A1","DUSP1","TAGLN2","FOS","TAGAP","S100A10")

quantile(Wbcross[,1],probs = seq(0,1,0.01))
Wbcross[intersect(keyrepgenes,rownames(Wbcross)),]

#######################
#### USAE: Doing the same but unsupervised

latent_size = 50

enc_input = layer_input(shape = ncol(X.train))
enc_output = enc_input %>% 
  layer_gaussian_noise(0.01) %>%
  layer_dropout(rate=0.5) %>%
  layer_dense(units=150, kernel_constraint = constraint_maxnorm(3)) %>% 
  layer_activation_leaky_relu(alpha=0.1) %>%
  layer_dense(units=latent_size,kernel_constraint = constraint_maxnorm(3)) %>%
  layer_activation_leaky_relu(alpha=0.1)

encoder = keras_model(enc_input, enc_output,name="encoder")

dec_input = layer_input(shape = latent_size)
dec_output = dec_input %>% 
  layer_dense(units=150) %>% 
  layer_activation_leaky_relu(alpha=0.1) %>% 
  layer_dense(units = ncol(X.train), name = "reconstruction_output") %>% 
  layer_activation_leaky_relu(alpha=0.1) 

decoder = keras_model(dec_input, dec_output,name="decoder")
ann_input = layer_input(shape = latent_size)
auxiliary_output <- ann_input %>% 
  layer_dense(units=latent_size) %>%
  layer_activation_leaky_relu(0.1) %>%
  layer_dense(units=ncol(Y.train),activation="linear")
ann <- keras_model(ann_input, auxiliary_output,name="ANN")

aen_input = layer_input(shape = ncol(X.train))
aen_output1 = aen_input %>% 
  encoder() %>% 
  decoder()

aen_output2 = aen_input %>% 
  encoder() %>% 
  ann()

aen = keras_model(aen_input, c(aen_output1,aen_output2),name="Modelfinal")
lr_schedule = optimizer_adam(
  learning_rate=0.001, ### very small
  beta_1 =0.9)
aen %>% compile(
  optimizer = lr_schedule,
  loss = 'mse',
  loss_weights = c(1,0), # <- meant to run unsupervised AE
  metrics = c("mse")
)
history <- aen %>% fit(X.train,list(X.train,Y.train), epochs=30, batch_size=64, validation_split=0.2) 


# 30 steps were sufficient!
history$params

#train data 
latentspace.train.us = encoder %>% predict(X.train)
decoded_space = decoder %>% predict(latentspace.train.us)
decoded_space[1:5,1:4]
plot(decoded_space[,1],X.train[,1]);abline(a=0,b=1)
plot(decoded_space[,10],X.train[,10]);abline(a=0,b=1)

BPpred.train.us <- ann %>% predict(latentspace.train.us)
plot(BPpred.train.us[,1],Y.train[,1],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
plot(BPpred.train.us[,2],Y.train[,2],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
cor(BPpred.train.us[,1],Y.train[,1])
cor(BPpred.train.us[,2],Y.train[,2])



# test data
latentspace.test.us = encoder %>% predict(X.test)
decoded_space = decoder %>% predict(latentspace.test.us)
latentspace.test.us[1:5,1:4]
plot(decoded_space[,1],X.test[,1]);abline(a=0,b=1)
plot(decoded_space[,10],X.test[,10]);abline(a=0,b=1)

BPpred.test.us <- ann %>% predict(latentspace.test.us)
plot(BPpred.test.us[,1],Y.test[,1],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
plot(BPpred.test.us[,2],Y.test[,2],xlim=c(-3,3),ylim=c(-3,3));abline(a=0,b=1)
cor(BPpred.test.us[,1],Y.test[,1])
cor(BPpred.test.us[,2],Y.test[,2])
summary(lm(BPpred.test.us[,1]~Y.test[,1]-1))
summary(lm(BPpred.test.us[,2]~Y.test[,2]-1))


### checking spaces before saving

dim(latentspace.test.us)
dim(latentspace.test.ss)
dim(latentspace.train.us)
dim(latentspace.train.ss)

save(latentspace.train.ss, latentspace.train.us,
     latentspace.test.ss,latentspace.test.us,
     file = "SubspaceTranscriptomics.RData")
