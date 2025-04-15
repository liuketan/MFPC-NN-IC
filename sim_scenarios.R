library(survival)
library(icenReg)
library(ICcforest)
library(dplyr)
library(ggplot2)
library(reshape2)
library(caret)
library(intcensROC)
library(keras)
library(tensorflow)
library(funData)
library(MFPCA)
library(MASS)
library(plyr)

source("fun_NN-IC.R")
source("MFPCA_fun.r")
source("data_gen.r")
source("cal_AUC_BS.r")



Scene = "scenario1"

epoch <- 500
batch_size <- 10
num_nodes <- 60
string_activation <- "selu"
num_l1 <- 0.01
num_dropout <- 0.5
num_lr <- 0.0002
num_layer <- 2

num_m = 3
num_l = 0
num_u = 14

#mtry = sqrt(5)
mtry = cores = NULL

I = 2000
I.test = 700
J = 8

nsim = 100
nbasis = 4
obstime = c(0,0.5,1,1.5,2,2.5,3,3.5)
Mu = 0.3 
argvals = obstime / max(obstime) # scale it to 0-1
pred_tstar = c(0.5,1)  # predict time t*= 0.5 and 1 
delta_t = c(0.5,1)     # predict time window delta_t = 0.5 and 1, meaning that we use the information at time 0.5 and 1 to predict 1, 1.5 and 2

k=2
lambda=1

Brier_matrix = Brier_f_matrix = Brier_icen_matrix = AUC_icen_matrix = AUC_matrix = AUC_f_matrix = True_AUC_matrix = array(NA,dim = c(nsim,length(pred_tstar),2),dimnames = list(NULL,c(paste0("tstar=",pred_tstar)),c(paste0("delta_t=",delta_t))))
BS = AUC = BS_sd = AUC_sd = array(NA,dim = c(3,length(pred_tstar),2),dimnames = list(c("NN","ICcforest","IcenReg"),c(paste0("tstar=",pred_tstar)),c(paste0("delta_t=",delta_t))))
TRUE_AUC = matrix(NA,nrow = 2,ncol = length(pred_tstar))
rownames(TRUE_AUC) = paste0("delta_t=",delta_t)
Censor_rate = Interval_length = c()
set.seed(888)

for (irun in 1:nsim) {
  
  tryCatch({
    print(c("###############################################################################################",
            paste("###########------------------------- This is the",irun,"th ------------------------------###########",sep = " "),
            "###############################################################################################"))
    
    if (Scene == "scenario1"){
      data = sim_multi_linear(I,J,MU = Mu, obstime, miss = F)
    }else if (Scene == "scenario2"){
      data = sim_multi_nonlinear(I,J,MU = Mu, obstime, miss = F)
    }else if(Scene == "scenario3"){
      data = sim_ht_nonlinear(I,J,MU = Mu, obstime, miss = F)
    }else if(Scene == "scenario4"){
      data = sim_both_nonlinear(I,J,MU = Mu, obstime, miss = F)
    }
    
    long = data$long
    surv = data$surv  
    min(surv$time)
    Interval_length[irun] = mean(surv$Right[!is.infinite(surv$Right)]-surv$Left[!is.infinite(surv$Right)])
    Censor_rate[irun] = sum(surv$event==0)/I 
    long_array = array(NA, c(nrow(surv), length(obstime), 5))
    
    visits = c()
    for (i in 1:nrow(surv)) {
      visits[i] = sum(long$id==sort(surv$id)[i])
    }
    for(i in 1:nrow(surv)){
      long_array[i,1:visits[i], 1] = long$Y1[long$id == surv$id[i]]
      long_array[i,1:visits[i], 2] = long$Y2[long$id == surv$id[i]]
      long_array[i,1:visits[i], 3]= long$Y3[long$id == surv$id[i]]
      long_array[i,1:visits[i], 4] = long$Y4[long$id == surv$id[i]]
      long_array[i,1:visits[i], 5]= long$Y5[long$id == surv$id[i]]
    }
    
    test_idx <- sample(1:I,size = I.test,replace = FALSE)
    surv.test = surv[test_idx,]
    surv.train = surv[-test_idx,]
    
    ## scale only use the training sd
    sd = c()
    long_scale = long_array
    for (z in 1:5){
      sd[z] = sd(long_array[test_idx,,z],na.rm = TRUE)
      long_scale[,,z] = long_array[,,z] / sd[z]
      print(var(as.vector(long_scale[,,z]),na.rm = TRUE)) 
    }
    
    long.test = long_scale[test_idx,,]
    long.train = long_scale[-test_idx,,]
    
    ############ MFPCA ###########
    
    test_idx <- sample(1:I,size = I.test,replace = FALSE)
    surv.test = surv[test_idx,]
    surv.train = surv[-test_idx,]
    
    
    # univariate FPCA via PACE
    Xi.train = L = phi.train = meanFun.train =  NULL
    
    for(p in 1:5){
      tmp.ufpca = uPACE(long.train[,,p], argvals, nbasis=nbasis)
      Xi.train = cbind(Xi.train, tmp.ufpca$scores) # FPC scores
      L = c(L, dim(tmp.ufpca$scores)[2])
      phi.train[[p]] = t(tmp.ufpca$functions@X) # FPC eigenfunctions
      meanFun.train[[p]] = tmp.ufpca$mu@X # estimated mean functions
    }
    
    # multivariate FPCA
    mFPCA.train = mFPCA(Xi=Xi.train, phi=phi.train, p=5, L=L,nsurv = nrow(surv.train) )
    rho.train = mFPCA.train$rho[,1:nbasis]  #MFPC scores
    pve = mFPCA.train$pve
    psi = mFPCA.train$psi
    Cms = mFPCA.train$Cms
    colnames(rho.train) = paste0("rho",1:nbasis)
    surv.all_train = data.frame(surv.train,rho.train)
    
    left_dat = as.matrix(surv.all_train$Left)
    right_dat = as.matrix(surv.all_train$Right)
    status_dat = c()
    for (z in 1:nrow(left_dat)) {
      status_dat[z] = ifelse(is.infinite(right_dat[z]),0,1)
    }
    status_dat = as.matrix(status_dat)
    pred = surv.all_train[,7:ncol(surv.all_train)]
    pred = as.matrix(pred)
    right_dat[is.infinite(right_dat)] = max(surv.all_train$Right[is.finite(surv.all_train$Right)]) + 0.1
    right_dat = as.matrix(right_dat)
    rm(model)
    model <- build_model_ic(left_dat, right_dat, pred, m = num_m, l = num_l, u = num_u, num_nodes, string_activation, num_l1, num_dropout, num_lr, num_layer)
    
    model %>% fit(
      list(left_dat, right_dat, pred),
      status_dat,
      epochs = epoch,
      batch_size = batch_size,
      verbose = 1
    )
    
    surv.fore_train <- surv.all_train
    surv.fore_train$Right <- ifelse(is.infinite(surv.fore_train$Right),999999,surv.fore_train$Right)
    model2 <- ICcforest::ICcforest(formula = Surv(Left,Right, type="interval2") ~ W1+W2+W3+rho1+rho2+rho3+rho4, 
                                   data = surv.fore_train, cores = cores, mtry = mtry)
    
    model_icen <- ic_sp(Surv(Left, Right, type = "interval2") ~ W1+W2+W3+rho1+rho2+rho3, model = "ph", data = surv.all_train)
    
    ith = 0
    for(t in pred_tstar){
      ith = ith+1
      lp1 = lp2 = lp3 = c()
      ##delete the t=0, so first t = 0.5
      tmp.id = surv.test[surv.test$time>t, "id"]  # subjects still event-free at landmark time 
      tmp.surv.data = surv.test[surv.test$time>t, ] # filter the data
      tmp.data = long_array[tmp.id, , ] # subset longitudinal outcomes
      tmp.data[,(which(t==obstime)+1):length(obstime),] = NA  # retain the measurements before landmark time
      
      # univariate FPC 
      Xi.test = NULL
      for(p in 1:5){
        tmp.ufpca = uPACE(long.train[,,p], argvals, tmp.data[,,p], nbasis=nbasis)
        Xi.test = cbind(Xi.test, tmp.ufpca$scores) # dynamic FPC scores for test subjects 
      }
      
      # estimate MFPC scores for test subjects
      rho.test = mfpca.score(Xi.test, Cms)[,1:nbasis]
      colnames(rho.test) =  paste0("rho",1:nbasis)
      
      # prepare test data 
      surv.all_test = data.frame(tmp.surv.data,rho.test)
      surv.fore_test <- surv.all_test # ICcforest
      surv.fore_test$Right <- ifelse(is.infinite(surv.fore_test$Right),999999,surv.fore_test$Right)
      left_test = as.matrix(surv.all_test$Left) # NN-IC
      right_test = as.matrix(surv.all_test$Right)
      status_test = c()
      for (z in 1:nrow(left_test)) {
        status_test[z] = ifelse(is.infinite(right_test[z]),0,1)
      }
      right_test[is.infinite(right_test)] = max(surv.all_test$Right[is.finite(surv.all_test$Right)]) + 0.1
      pred_test <- surv.all_test[,7:ncol(surv.all_test)]
      pred_test <- as.matrix(pred_test)
      
      # linear predictor to calculate AUC
      lp1 <- predict(model_icen,newdata = surv.all_test,type = "lp")
      
      Pred <- predict(model2, type = "prob", newdata = surv.fore_test, OOB = F)
      for (z in 1:nrow(surv.fore_test)) {
        lp2[z] = Pred[[z]]$llk
      }
      
      lp3 <- model %>% predict(list(left_test, right_test, pred_test))
      lp3 <- lp3[,3] # outputpred，g(z)
      
      lp_t <- as.numeric(log(-log(surv.all_test$true.prob[,ith])/ exp(-7)))
      
      for (d in 1:length(delta_t)) {
        dt = delta_t[d]
        ### MFPC-icenReg ####
        Brier_icen_matrix[irun,ith,d] <- sbrier_IC_updated(obj = survival::Surv(surv.fore_test$Left,surv.fore_test$Right, type = "interval2"),
                                                           pred = model_icen,
                                                           newData = surv.fore_test, # for icenReg, need to specify newdata
                                                           tstar = t,
                                                           tp = t+dt)
        AUC_icen_matrix[irun,ith,d] <- cal.survAUC(lp1,surv.all_test,t+dt)
        
        ### MFPC-ICcforest ###
        Brier_f_matrix[irun,ith,d] <- sbrier_IC_updated(obj = survival::Surv(surv.fore_test$Left, surv.fore_test$Right, type = "interval2"),
                                                        pred = Pred, 
                                                        tstar = t,
                                                        tp = t+dt)
        AUC_f_matrix[irun,ith,d] <- cal.survAUC(lp2,surv.all_test,t+dt)
        
        ### MFPC-NN-IC ###
        brier<-  brier_ic_nn(obj = model,
                             left_new = left_test,
                             right_new = right_test,
                             x_new = pred_test,
                             tstar = t,
                             t = t+dt,
                             type = "BS")
        Brier_matrix[irun,ith,d] <- brier[[length(brier)]]
        AUC_matrix[irun,ith,d] <- cal.survAUC(lp3,surv.all_test,t+dt)
        
        ##### true AUC ####
        True_AUC_matrix[irun,ith,d] = cal.survAUC(lp_t,surv.all_test,t+dt)
      }
    }
    
  },
  error=function(e){cat("Error",conditionMessage(e), "\n")})
}


calculate_stats <- function(matrix_list, stat_fun) {
  result <- array(NA, dim = c(3, dim(matrix_list[[1]])[2], 2)) # 3 methods, ncol, 2 time points
  for (i in 1:3) {
    for (j in 1:2) {
      result[i,,j] <- apply(na.omit(matrix_list[[i]][,,j]), 2, stat_fun)
    }
  }
  return(result)
}


BS_matrices <- list(Brier_matrix, Brier_f_matrix, Brier_icen_matrix)
AUC_matrices <- list(AUC_matrix, AUC_f_matrix, AUC_icen_matrix)
TRUE_AUC <- sapply(1:2, function(j) apply(na.omit(True_AUC_matrix[,,j]), 2, mean))

rm(sd)
BS <- calculate_stats(BS_matrices, mean)
AUC <- calculate_stats(AUC_matrices, mean)
BS_sd <- calculate_stats(BS_matrices, sd)
AUC_sd <- calculate_stats(AUC_matrices, sd)

output_matrix <- matrix(NA, nrow = length(pred_tstar) * length(delta_t), ncol = 9)
output_matrix[,1] <- rep(pred_tstar, each = 2)
output_matrix[,2] <- rep(delta_t, 2)
output_matrix[,3] <- round(as.vector(TRUE_AUC), 3)

format_stats <- function(AUC, AUC_sd, BS, BS_sd) {
  paste0(round(as.vector(t(AUC)), 3), " (", round(as.vector(AUC_sd), 3), ")")
}

output_matrix[,4] <- format_stats(AUC[1,,], AUC_sd[1,,], BS[1,,], BS_sd[1,,])
output_matrix[,5] <- format_stats(BS[1,,], BS_sd[1,,], AUC[1,,], AUC_sd[1,,])
output_matrix[,6] <- format_stats(AUC[2,,], AUC_sd[2,,], BS[2,,], BS_sd[2,,])
output_matrix[,7] <- format_stats(BS[2,,], BS_sd[2,,], AUC[2,,], AUC_sd[2,,])
output_matrix[,8] <- format_stats(AUC[3,,], AUC_sd[3,,], BS[3,,], BS_sd[3,,])
output_matrix[,9] <- format_stats(BS[3,,], BS_sd[3,,], AUC[3,,], AUC_sd[3,,])

colnames(output_matrix) <- c("t", "dt", "True AUC", "NN-IC AUC", "NN-IC BS", "ICcForest AUC", "ICcForest BS", "icenReg AUC", "icenreg BS")
write.csv(output_matrix, file = "Sim_1.csv")


