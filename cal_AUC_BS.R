##################### AUC BS calculation functions ######################

# brier_ic_nn        # calculate dynamic Brier Score for NN-IC
# sbrier_IC_updated  # calculate Brier Score for interval-censored data (icenReg and ICcforest)
# cal.survAUC        # calculate AUC for interval-censored data

brier_ic_nn <- function(obj, left_new, right_new, x_new, tstar,tp,btime = NULL, type = c("IBS", "BS")) {
  
  N <- nrow(left_new)
  intL <- left_new
  intR <- right_new
  times <- c(intL, intR)
  if (is.null(btime)) btime=range(tstar, tp)
  
  # the observation time in test set
  newtime <- unique(times)
  newtime <- newtime[order(newtime)]
  
  ltstar <- max(which(newtime <= tstar))
  utstar <- which(newtime >= tstar)[1]
  ltstar <- newtime[ltstar]
  utstar <- newtime[utstar]
  truetstar <- NULL
  if(abs(ltstar - tstar) < abs(utstar - tstar)){
    truetstar <- ltstar
  }else{truetstar <- utstar}
  truetstaridx = which(times == truetstar)[1]
  
  ltp <- max(which(newtime <= tp))
  utp <- which(newtime >= tp)[1]
  ltp <- newtime[ltp]
  utp <- newtime[utp]
  truetp <- NULL
  if(abs(ltp - tp) < abs(utp - tp)){
    truetp <- ltp
  }else{truetp <- utp}
  truetpidx = which(times == truetp)[1]
  
  # Calculate S via NN-IC output
  out.1 <- obj %>% predict(list(left_new, right_new, x_new))
  out.2 <- obj %>% predict(list(times, times, rbind(x_new, x_new)))
  S_times = t(k_get_value(exp(-1 * (k_dot(k_cast(matrix(exp(out.1[,3]), ncol = 1), dtype = tf$float32), k_cast(matrix(out.2[,1], nrow = 1), dtype = tf$float32))))))
  
  bsc_temp <- matrix(0, nrow = N, ncol = 1)
  k = 0
  sumIt = 0
  for (j in 1:N){
    # Estimate It = (St0-Sr)/(Sl-Sr)
    St0 <- S_times[truetstaridx,j]
    if (times[truetstaridx] <= intL[j]){
      It = 1
    } else if(times[truetstaridx] > intR[j]){
      It = 0
    } else {
      Sr <- S_times[N + j, j]
      Sl <- S_times[j, j]
      if (Sl != Sr){
        It = (St0 - Sr)/(Sl - Sr)
      } else {
        It = NULL
      }
    }
    if(!is.null(It)){sumIt = as.numeric(sumIt) + as.numeric(It)}
    
    # Estimate IY = (St-Sr)/(Sl-SR)
    Stp <- S_times[truetpidx,j]
    St <- Stp/St0
    if (times[truetpidx] <= intL[j]){
      IY = 1
    } else if(times[truetpidx] > intR[j]){
      IY = 0
    } else {
      Sr <- S_times[N + j, j]
      Sl <- S_times[j, j]
      if (Sl != Sr){
        IY = (Stp - Sr)/(Sl - Sr)
      } else {
        IY = NULL
      }
    }
    if (is.null(IY)==0){
      k = k + 1
      bsc_temp[k] = It * (IY - St)^2
    }
  }
  
  bsc = sum(bsc_temp[1:k]) / sumIt
  bsc
  
}


sbrier_IC_updated <- function(obj, pred, tstar, tp, newData = NULL) {
  # obj = survival::Surv(data$Left, data$Right, type = "interval2")
  # pred = model1
  # tstar = 0.5
  # tp = 1
  # newData = data_2 # for icenReg, need to specify newdata
  
  
  
  if (!inherits(obj, "Surv"))
    stop("obj is not of class Surv")
  class(obj) <- NULL
  # number of obs
  N <- nrow(obj)
  
  intL <- obj[, 1]
  intR <- obj[, 2]
  # the observation time in test set
  time <- c(intL,intR)
  newtime <- unique(time)
  newtime <- newtime[order(newtime)]
  
  ltstar <- max(which(newtime <= tstar))
  utstar <- which(newtime >= tstar)[1]
  ltstar <- newtime[ltstar]
  utstar <- newtime[utstar]
  truetstar <- NULL
  if(abs(ltstar - tstar) < abs(utstar - tstar)){
    truetstar <- ltstar
  }else{truetstar <- utstar}
  truetstaridx = which(time == truetstar)[1]
  
  ltp <- max(which(newtime <= tp))
  utp <- which(newtime >= tp)[1]
  ltp <- newtime[ltp]
  utp <- newtime[utp]
  truetp <- NULL
  if(abs(ltp - tp) < abs(utp - tp)){
    truetp <- ltp
  }else{truetp <- utp}
  truetpidx = which(time == truetp)[1]
  
  
  ptype <- class(pred)[[1]] # ic_np, ic_ph, ic_po
  
  
  survs <- NULL
  switch(ptype, survfit = {
    survs <- try(ipred::getsurv(pred, time), silent=T)
    if(class(survs) == "try-error"){
      stop("please check whether the package 'ipred' is installed")
    }
    survs <- matrix(rep(survs, N), nrow = length(time))
  }, ic_np = {
    
    if (is.null(newData)) {
      survs <- try(getFitEsts_Surv(Curve = pred, teval = time), silent=T) # newData = newData
      if(class(survs) == "try-error"){
        stop("please check whether the package 'icenReg' is installed")
      }
      survs <- matrix(rep(survs, N), nrow = length(time))
    } else {
      survs <- lapply(time, function(x) {getFitEsts_Surv(Curve = pred, teval = x, newData = newData)})
      survs <- t(matrix(unlist(survs), ncol = length(time), nrow = nrow(newData)))
    }
    
  }, ic_ph = {
    
    if (is.null(newData)) {
      survs <- try(getFitEsts_Surv(Curve = pred, teval = time), silent=T) # newData = newData
      if(class(survs) == "try-error"){
        stop("please check whether the package 'icenReg' is installed")
      }
      survs <- matrix(rep(survs, N), nrow = length(time))
    } else {
      survs <- lapply(time, function(x) {getFitEsts_Surv(Curve = pred, teval = x, newData = newData)})
      survs <- t(matrix(unlist(survs), ncol = length(time), nrow = nrow(newData)))
    }
    
  }, ic_po = {
    
    if (is.null(newData)) {
      survs <- try(getFitEsts_Surv(Curve = pred, teval = time), silent=T) # newData = newData
      if(class(survs) == "try-error"){
        stop("please check whether the package 'icenReg' is installed")
      }
      survs <- matrix(rep(survs, N), nrow = length(time))
    } else {
      survs <- lapply(time, function(x) {getFitEsts_Surv(Curve = pred, teval = x, newData = newData)})
      survs <- t(matrix(unlist(survs), ncol = length(time), nrow = nrow(newData)))
    }
    
  },list = {
    if (!inherits(pred[[1]], "survfit") && !inherits(pred[[1]], "ic_ph") && !inherits(pred[[1]], "ic_po") && !inherits(pred[[1]], "ic_np"))
      stop("pred is not a list of survfit/ic_ph/ic_np/ic_po objects; \n  if pred is an ICcforest prediction object, please make sure it is created using the setting type = 'prob'")
    if (length(pred) != N) stop("pred must be of length(time)")
    if (inherits(pred[[1]], "survfit")){
      M.list <- try(lapply(pred, ipred::getsurv, times = time),silent=T)
      if (class(M.list) == "try-error"){
        stop("please check whether the package 'ipred' is installed")
      }
      survs <- matrix(unlist(lapply(M.list, function(x) x)),
                      nrow = length(time), ncol = N) ## each obs in one column
    } else if (inherits(pred[[1]], "ic_np") || (inherits(pred[[1]], "ic_ph")) || (inherits(pred[[1]], "ic_po"))) {
      M.list <- try(lapply(pred, getFitEsts_Surv, teval = time), silent=T) # no need for newData here since it is included in pred
      if(class(M.list) == "try-error"){
        stop("please check whether the package 'icenReg' is installed")
      }
      survs <- matrix(unlist(lapply(M.list, function(x) x)),
                      nrow = length(time), ncol = N) ## each obs in one column
    }
  }, vector = {
    if (length(pred) != N) stop("pred must be of length(time)")
    if (length(time) != 1) stop("cannot compute integrated Brier score with pred; \n  if pred is an ICcforest prediction object, please make sure it is created using the setting type = 'prob'")
    survs <- pred
  }, matrix = {
    if (all(dim(pred) == c(length(time), N))) survs <- pred
    else stop("wrong dimensions of pred; \n  if pred is an ICcforest prediction object, please make sure it is created using the setting type = 'prob'")
  })
  if (is.null(survs))
    stop("unknown type of pred")
  
  bsc_temp <- matrix(0, nrow = N, ncol = 1)
  k = 0
  sumIt = 0
  for (j in 1:N){
    # Estimate It = (St0-Sr)/(Sl-Sr)
    St0 <- survs[truetstaridx,j]
    if (time[truetstaridx] <= intL[j]){
      It = 1
    } else if(time[truetstaridx] > intR[j]){
      It = 0
    } else {
      Sr <- survs[N+j, j]
      Sl <- survs[j, j]
      if (Sl != Sr){
        It = (St0 - Sr)/(Sl - Sr)
      } else {
        It = NULL
      }
    }
    
    if(!is.null(It)){
      sumIt = as.numeric(sumIt) + as.numeric(It)
      
      # Estimate IY = (St-Sr)/(Sl-SR)
      Stp <- survs[truetpidx,j]
      St <- Stp/St0
      if (time[truetpidx] <= intL[j]){
        IY = 1
      } else if(time[truetpidx] > intR[j]){
        IY = 0
      } else {
        Sr <- survs[N+j, j]
        Sl <- survs[j, j]
        if (Sl != Sr){
          IY = (Stp - Sr)/(Sl - Sr)
        } else {
          IY = NULL
        }
      }
      if (is.null(IY)==0){
        k = k + 1
        bsc_temp[k] = It * (IY - St)^2
      }
    }
  }
  bsc = sum(bsc_temp[!is.nan(bsc_temp)]) / sumIt
  bsc
}

### for icenreg,lp <- predict(model,data = surv.test,type = "lp")
### for iccforest:   
# for (z in 1:nrow(surv.test)) {
# lp[z] = Pred[[z]]$llk
# }
### for NN,   lp <- model %>% predict(list(left_test, right_test, pred_test)),lp = lp[,3]

cal.survAUC <- function(lp,surv_data,ptime){
  Marker = (lp - min(lp[!is.infinite(lp)]))/(max(lp) - min(lp[!is.infinite(lp)])) + 0.1 # need to be positive
  Marker = ifelse(Marker<=0,0.1,Marker)
  U <- as.vector(surv_data$Left)
  V <- as.vector(surv_data$Right)
  #U <- ifelse(surv_data$Left == 0, V, U) # for left-censor
  #Delta <- ifelse(surv_data$Left == 0, 1, 2) # for left-censor
  V <- ifelse(is.infinite(surv_data$Right), U, V) # for right-censor
  Delta <- ifelse(is.infinite(surv_data$Right), 3, 2) # 3 for right-censor, 2 for interval-censor
  res <- intcensROC(U, V, Marker, Delta, ptime, gridNumber = 500)
  AUCi <- intcensAUC(res)
  return(AUCi)
}
