################################## Model Compare ############################################
########################### ADNI MCI + CN 1702 individuals ##################################
### Model 1. Baseline scores NN-IC: 4 baseline + 5 baseline cognitive scores
### Model 2. MFPC-NN-IC: 4 baseline + 5 tvc 
### Model 3. MFPC-icenReg: 4 baseline + 5 tvc 
### Model 4. MFPC-ICcforest: 4 baseline + 5 tvc 
#############################################################################################

library(reticulate)
library(tensorflow)
library(keras)
library(MASS)
library(icenReg)
library(ICcforest)
library(dplyr)
library(reshape2)
library(caret)
library(intcensROC)
library(survival)
library(funData)
library(MFPCA)
library(plyr)


tstar = c(0.5,1,1.5,2) # evaluate time point
delta_t = c(0.5,1)


## hyperparameter for NN-IC
k=5
epoch <- 1000
batch_size <- 60
num_nodes <- 45
string_activation <- "selu"
num_l1 <- 0.02
num_dropout <- 0.7
num_lr <- 0.0002
num_layer <- 1
num_m = 3
num_l = 0
num_u = 14

set.seed(66)

mtry = sqrt(5)
cores = NULL

nbasis = 3

source("functions.R")
source("fun_DNN-IC.R")
source("AUC_BS_fun.r")
data <- read.csv("ADNI_all_noAD(MCI+CN).csv")
rawdata <- read.csv("ADNIMERGE.csv")

data$PTGENDER[data$PTGENDER == "Female"] <- 1
data$PTGENDER[data$PTGENDER == "Male"] <- 0
data$PTGENDER <- as.numeric(data$PTGENDER)
minmax <- function(x){
  newx <- (x - min(x))/(max(x) - min(x))
}

#### Data Prepare
data_bl <- data.frame(data$RID,data$status,data$Left,data$Right,data$AGE,data$PTGENDER,data$PTEDUCAT,data$APOE4,data$ADAS13_bl,
                      data$RAVLT_immediate_bl,data$RAVLT_learning_bl,
                      data$FAQ_bl,data$MMSE_bl)
colnames(data_bl) <- c("RID","status","Left","Right","AGE","GENDER","EDUCATE","APOE4","ADAS13","RAVLT_imm","RAVLT_learn","FAQ","MMSE")
data_bl <- data_bl[order(data_bl$RID),]
data_bl <- apply(data_bl, 2, as.numeric)
data1 <- cbind(data_bl[,1:4],apply(data_bl[,5:13], 2, minmax)) # for model2

rawdata <- rawdata[order(rawdata$RID),]
data_surv <- as.data.frame(data1[,1:8])

surv_dat <- data.frame(data_surv$RID,data_surv$Left,data_surv$Right,data_surv$AGE,data_surv$GENDER,data_surv$EDUCAT,data_surv$APOE4)
long_dat <- data.frame(rawdata$RID,rawdata$Years_bl,rawdata$VISCODE,rawdata$ADAS13,rawdata$MMSE,
                       rawdata$FAQ,rawdata$RAVLT_immediate,rawdata$RAVLT_learning)

for (i in 1:nrow(rawdata)) {
  for (j in 1:nrow(data_surv)) {
    if(rawdata$RID[i] == data_surv$RID[j]){
      rawdata$Left[i] = data_surv$Left[j]
      rawdata$Right[i] = data_surv$Right[j]
    }
  }
}

long_dat <- data.frame(rawdata$RID,rawdata$Years_bl,rawdata$VISCODE,rawdata$Right,rawdata$ADAS13,rawdata$MMSE,
                       rawdata$FAQ,rawdata$RAVLT_immediate,rawdata$RAVLT_learning)
surv_scale <- apply(surv_dat[4:7], 2, minmax)
surv_scale_dat <- cbind(surv_dat[1:3],surv_scale)

#### get 0.5 year gap ######
long_dat$rawdata.VISCODE = ifelse(long_dat$rawdata.VISCODE=="bl","m0",long_dat$rawdata.VISCODE)
long_dat$vis_month <- as.numeric(substring(long_dat$rawdata.VISCODE,2))
long_dat$vis_seq <- long_dat$vis_month/6+1
df_long <- subset(long_dat,!is.integer(vis_seq))
df_long <- subset(df_long,row.names(df_long)!=14117)
raw_long <- array(NA, c(nrow(data_surv), 32, 5)) # subject, vis_month, covariates
for (i in 1:nrow(data_surv)) {
  tmp_long = df_long[df_long$rawdata.RID==data_surv$RID[i],5:9]
  tmp_long$vis_seq = df_long[df_long$rawdata.RID==data_surv$RID[i],11]
  for (v in 1:length(tmp_long$vis_seq)) {
    raw_long[i,tmp_long$vis_seq[v],] = unlist(tmp_long[which(tmp_long$vis_seq==tmp_long$vis_seq[v]),1:5])
  }
}


Model1_AUC = Model2_AUC = Model3_AUC = Model4_AUC = 
  Model1_BS = Model2_BS = Model3_BS = Model4_BS = 
  array(NA, dim = c(5,length(tstar),2))

CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]  
  temp <- sample(n,datasize) 
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x]) 
  return(cvlist)
}

lp1 = lp2 = lp3 = lp4 = c()

tryCatch({
  for (a in 1:k){
    
    ################### data preparation ####################
    cvlist_ex <- CVgroup(k=5,datasize = nrow(data1)) # 5 fold
    data1_ex <- data.frame(data1[cvlist_ex[[a]],]) # 1 fold for external validation
    data1_in <- data1[-cvlist_ex[[a]],] # 4 fold for internal validation
    cvlist_in <- CVgroup(k=5,datasize = nrow(data1_in)) # internal 5 fold 
    data1_train <- data.frame(data1_in[-cvlist_in[[a]],] )# 4 folds for training
    base_data_ex = data1_ex[,1:8]
    
    #train
    left_dat = as.matrix(data1_train[, 3], nrow = nrow(data1_train))
    right_dat = as.matrix(data1_train[, 4], nrow = nrow(data1_train))
    status_dat = as.matrix(data1_train[, 2], nrow = nrow(data1_train))
    right_dat[is.infinite(right_dat)] = max(data1_train$Right[is.finite(data1_train$Right)]) + 0.1
    
    #ex
    left_ex = as.matrix(data1_ex[, 3], nrow = nrow(data1_ex))
    right_ex = as.matrix(data1_ex[, 4], nrow = nrow(data1_ex))
    status_ex = as.matrix(data1_ex[, 2], nrow = nrow(data1_ex))
    right_ex[is.infinite(right_ex)] = max(data1_ex$Right[is.finite(data1_ex$Right)]) + 0.1
    
    pred_train <- data1_train[,5:13]
    pred_train <- as.matrix(pred_train,ncol=9,nrow=nrow(data1_train))
    
    pred_ex <- data1_ex[,5:13]
    pred_ex <- as.matrix(pred_ex,ncol=9,nrow=nrow(data1_ex))
    
    ################### MFPCA data #######################
    longvar.ex <- raw_long[cvlist_ex[[a]],,]
    longvar.in <- raw_long[-cvlist_ex[[a]], , ]
    surv.in <- data_surv[-cvlist_ex[[a]],3:8]
    surv.ex <- data_surv[cvlist_ex[[a]],3:8]
    longvar.train = longvar.in[-cvlist_in[[a]],,]
    
    ### scale only use the training sd
    sd = c()
    long_scale = raw_long
    for (z in 1:5){
      sd[z] = sd(longvar.train[,,z],na.rm = TRUE)
      long_scale[,,z] = raw_long[,,z] / sd[z]
      print(var(as.vector(long_scale[,,z]),na.rm = TRUE)) 
    }
    
    long_scale_in = long_scale[-cvlist_ex[[a]],,]
    long_scale_ex = long_scale[cvlist_ex[[a]],,]
    longscale.train = long_scale_in[-cvlist_in[[a]],,]
    longscale.intest = long_scale_in[cvlist_in[[a]],,]
    surv.train = surv.in[-cvlist_in[[a]],]
    surv_in_test = surv.in[cvlist_in[[a]],]
    
    colnames(surv.ex) = colnames(surv.train) =  c("Left","Right","age","gender","educate","apoe4")
    
    argvals = seq(from = 0,to = 0.5*31,by = 0.5)
    Xi.train = L = phi.train = meanFun.train =  NULL
    
    ### Conduct MFPCA in training sets
    for(p in 1:5){
      tmp.ufpca = uPACE(longscale.train[,,p], argvals, nbasis=nbasis)
      Xi.train = cbind(Xi.train, tmp.ufpca$scores) # FPC scores
      L = c(L, dim(tmp.ufpca$scores)[2])
      phi.train[[p]] = t(tmp.ufpca$functions@X) # FPC eigenfunctions
      meanFun.train[[p]] = tmp.ufpca$mu@X # estimated mean functions
    }
    # MFPCA
    mFPCA.train = mFPCA(Xi=Xi.train, phi=phi.train, p=5, L=L,nsurv = nrow(surv.train) )
    rho.train = mFPCA.train$rho[,1:nbasis]  #MFPC scores
    pve = mFPCA.train$pve #propoetion of variance
    psi = mFPCA.train$psi
    Cms = mFPCA.train$Cms
    
    colnames(rho.train) = paste0("rho",1:3)
    surv.all_train = data.frame(surv.train,rho.train)
    
    #########################################################################
    ###########################  Training  ##################################
    #########################################################################
    
    ####### Model 1
    rm(Base_score_NN)
    Base_score_NN <- build_model_ic(left_dat, right_dat, pred_train, m = num_m, l = num_l, u = num_u, num_nodes, string_activation, num_l1, num_dropout, num_lr, num_layer)
    
    Base_score_NN %>% fit(
      list(left_dat, right_dat, pred_train),
      status_dat,
      epochs = epoch,
      batch_size = batch_size,
      verbose = 1
    )
    
    
    
    ###### Model 2 MFPC-NN-IC
    left_dat = as.matrix(surv.all_train$Left)
    right_dat = as.matrix(surv.all_train$Right)
    status_dat = c()
    for (z in 1:nrow(left_dat)) {
      status_dat[z] = ifelse(is.infinite(right_dat[z]),0,1)
    }
    status_dat = as.matrix(status_dat)
    pred = surv.all_train[3:9]
    pred = as.matrix(pred)
    right_dat[is.infinite(right_dat)] = max(surv.all_train$Right[is.finite(surv.all_train$Right)]) + 0.1
    right_dat = as.matrix(right_dat)
    rm(model_NN)
    model_NN <- build_model_ic(left_dat, right_dat, pred, m = num_m, l = num_l, u = num_u, num_nodes, string_activation, num_l1, num_dropout, num_lr, num_layer)
    
    model_NN %>% fit(
      list(left_dat, right_dat, pred),
      status_dat,
      epochs = epoch,
      batch_size = batch_size,
      verbose = 1
    )
    
    ###### Model 3 MFPC-Icenreg 
    model_icen <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ., model = "ph", data = surv.all_train)
    
    ###### Model 4 MFPC-ICcforest 
    fsurv.train = surv.all_train
    fsurv.train$Right <- ifelse(is.infinite(fsurv.train$Right),999,fsurv.train$Right)
    model_forest <- ICcforest::ICcforest(formula = survival::Surv(Left,Right, type="interval2") ~ ., 
                                         data = fsurv.train, cores = cores, mtry = mtry)
    
    
    ################# Test #####################
    for (i in 1:length(tstar)) {
      #### MFPCA test
      long_ex_NA = long_scale_ex
      # trajectory beyond tstar should be removed
      long_ex_NA[,(1+i*2):32,] = NA
      # univariate FPC 
      outer.test = NULL
      for(p in 1:5){
        outer.ufpca = uPACE(longscale.train[,,p], argvals, long_ex_NA[,,p], nbasis=nbasis)
        outer.test = cbind(outer.test, outer.ufpca$scores) # dynamic FPC scores for test subjects 
      }
      # estimate MFPC scores for test subjects
      rho.out = mfpca.score(outer.test, Cms)[,1:nbasis]
      colnames(rho.out) = paste0("rho",1:3)
      surv.ex_test = data.frame(surv.ex,rho.out)
      
      #test for MFPC-NN-IC
      left_test = as.matrix(surv.ex_test$Left)
      right_test = as.matrix(surv.ex_test$Right)
      status_ex_test = c()
      for (z in 1:nrow(left_test)) {
        status_ex_test[z] = ifelse(is.infinite(right_test[z]),0,1)
      }
      right_test[is.infinite(right_test)] = max(surv.ex_test$Right[is.finite(surv.ex_test$Right)]) + 0.1
      pred_test <- surv.ex_test[,3:9]
      pred_test <- as.matrix(pred_test)
      
      for (dt in 1:length(delta_t)) {
        final_t = tstar[i]+delta_t[dt]
        ######## Baseline scores NN-IC ##########
        brier<-  brier_ic_nn(obj = Base_score_NN,
                             left_new = left_ex,
                             right_new = right_ex,
                             x_new = pred_ex,
                             tstar = tstar[i],
                             t = final_t,
                             type = "BS")
        Model1_BS[a,i,dt] <- brier[[length(brier)]]
        
        lp1 <- Base_score_NN %>% predict(list(left_ex, right_ex, pred_ex))
        lp1 = lp1[,3]
        Model1_AUC[a,i,dt] = cal.survAUC(lp1,data1_ex,final_t)
        
        ######## MFPC-NN-IC ###########
        brier<-  brier_ic_nn(obj = model_NN,
                             left_new = left_ex,
                             right_new = right_ex,
                             x_new = pred_test,
                             tstar = tstar[i],
                             t = final_t,
                             type = "BS")
        Model2_BS[a,i,dt] <- brier[[length(brier)]]
        lp2 <- model_NN %>% predict(list(left_test, right_test, pred_test))
        lp2 = lp2[,3]
        Model2_AUC[a,i,dt] <- cal.survAUC(lp2,surv.ex_test,final_t)
        
        ########### MFPC-icenReg ##########
        Model3_BS[a,i,dt] = sbrier_IC_updated(obj = survival::Surv(surv.ex_test$Left, surv.ex_test$Right, type = "interval2"),
                                              pred = model_icen,
                                              newData = surv.ex_test,
                                              tstar = tstar[i],
                                              tp = final_t)
        lp3 <- predict(model_icen,newdata = surv.ex_test,type = "lp")
        Model3_AUC[a,i,dt] <- cal.survAUC(lp3,surv.ex_test,final_t)
        
        ########### MFPC-ICcforest ############
        fsurv.ex = surv.ex_test
        fsurv.ex$Right <- ifelse(is.infinite(fsurv.ex$Right),999,fsurv.ex$Right)
        Pred <- predict(model_forest, type = "prob", newdata = fsurv.ex, OOB = F)
        Model4_BS[a,i,dt] <- sbrier_IC_updated(obj = survival::Surv(fsurv.ex$Left, fsurv.ex$Right, type = "interval2"),
                                               pred = Pred, 
                                               tstar = tstar[i],
                                               tp = final_t)
        for (z in 1:nrow(fsurv.ex)) {lp4[z] = Pred[[z]]$llk}
        lp4 = as.numeric(lp4)
        Model4_AUC[a,i,dt] <- cal.survAUC(lp4,fsurv.ex,final_t)
      }
    }
  }
},
error=function(e){cat("Error",conditionMessage(e), "\n")})

AUC1 = AUC2 = AUC3 = AUC4 = BS1 = BS2 = BS3 = BS4 =  
  matrix(NA,nrow = length(tstar),ncol = length(delta_t),dimnames = list(c(),c("dt=0.5","dt=1")))
for(i in 1:2){
  AUC1[,i] = apply(Model1_AUC[,,i],2,mean)
  AUC2[,i] = apply(Model2_AUC[,,i],2,mean)
  AUC3[,i] = apply(Model3_AUC[,,i],2,mean)
  AUC4[,i] = apply(Model4_AUC[,,i],2,mean)
  BS1[,i] = apply(Model1_BS[,,i],2,mean)
  BS2[,i] = apply(Model2_BS[,,i],2,mean)
  BS3[,i] = apply(Model3_BS[,,i],2,mean)
  BS4[,i] = apply(Model4_BS[,,i],2,mean)
}

############################### Output AUC and BS ###############################
dt1_AUC = data.frame(tstar,c(rep(delta_t[1],length(tstar))),AUC1[,1],AUC2[,1],AUC3[,1],AUC4[,1])
dt2_AUC = data.frame(tstar,c(rep(delta_t[2],length(tstar))),AUC1[,2],AUC2[,2],AUC3[,2],AUC4[,2])
dt1_BS = data.frame(tstar,c(rep(delta_t[1],length(tstar))),BS1[,1],BS2[,1],BS3[,1],BS4[,1])
dt2_BS = data.frame(tstar,c(rep(delta_t[2],length(tstar))),BS1[,2],BS2[,2],BS3[,2],BS4[,2])
colnames(dt1_AUC) = colnames(dt2_AUC) = colnames(dt1_BS) = colnames(dt2_BS) = 
  c("tstar","dt","Baseline scores NN-IC","MFPC-NN-IC","MFPC icenReg","MFPC ICcforest")
AUC = BS = matrix(NA,nrow = length(tstar)*length(delta_t),ncol = 6)
AUC = rbind(dt1_AUC,dt2_AUC)
BS = rbind(dt1_BS,dt2_BS)
