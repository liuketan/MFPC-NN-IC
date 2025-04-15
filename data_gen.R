######################################################################################################
################################## Generate simulation data ##########################################
######################################################################################################

# sim_multi_linear(Scenario1)         # 3 fixed and 5 longitudinal, linear 
# sim_multi_nonlinear(Scenario2)      # only time-varying nonlinear
# sim_ht_nonlinear(Scenario3)         # only time-fixed nonlinear, interaction of W1 and W3 
# sim_both_nonlinear(Scenario4)       # both time-varying and time-fixed are nonlinear


sim_multi_linear = function(I, J, MU, obstime = obstime, miss = FALSE, miss.rate = 0.1){
  
  # I : number of subjects
  # J : number of visits
  # obstime: observation times
  # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
  # miss.rate: missing rate.
  
  N = I*J
  
  # longitudinal submodel  
  beta0 = c(1,1,1,2,1)
  beta1 = c(1.5,1,-1,1.5,1)
  betat = c(0.5,0.5,1.2,-1, 1.5) 
  ####generate subject-specific random effects from multi-normal distribution
  b.eta = runif(10,-0.5,0.5)
  b.sigma2 = c(0.5,1,1.5,2,2.5)
  e.var = c(1,1,1,1,1)
  bSigma = diag(b.sigma2)
  bSigma[1,2] = bSigma[2,1] = sqrt(b.sigma2[1]*b.sigma2[2])*b.eta[1]
  bSigma[1,3] = bSigma[3,1] = sqrt(b.sigma2[1]*b.sigma2[3])*b.eta[2]
  bSigma[1,4] = bSigma[4,1] = sqrt(b.sigma2[1]*b.sigma2[4])*b.eta[3]
  bSigma[1,5] = bSigma[5,1] = sqrt(b.sigma2[1]*b.sigma2[5])*b.eta[4]
  bSigma[2,3] = bSigma[3,2] = sqrt(b.sigma2[2]*b.sigma2[3])*b.eta[5]
  bSigma[2,4] = bSigma[4,2] = sqrt(b.sigma2[2]*b.sigma2[4])*b.eta[6]
  bSigma[2,5] = bSigma[5,2] = sqrt(b.sigma2[2]*b.sigma2[5])*b.eta[7]
  bSigma[3,4] = bSigma[4,3] = sqrt(b.sigma2[3]*b.sigma2[4])*b.eta[8]
  bSigma[3,5] = bSigma[5,3] = sqrt(b.sigma2[3]*b.sigma2[5])*b.eta[9]
  bSigma[4,5] = bSigma[5,4] = sqrt(b.sigma2[4]*b.sigma2[5])*b.eta[10]
  
  
  # sample covariate
  X = rep(rnorm(I, 6, 1), each=J)
  # sample random effect
  ranef = mvrnorm(I, c(0,0,0,0,0), bSigma)
  id = rep(1:I,each=J)
  ranef = ranef[id,]
  # construct longitudinal submodel
  eta.long = matrix(0, nrow=N, ncol=5)
  for(i in 1:5){
    eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
  }
  
  
  
  # survival submodel
  alpha = c(0.1,0.07,0.15,0.25,0.1)# longitudinal coefficients
  gamma <- c(2,-2,-1.5)
  W1 = W2 = W3 = W = NA
  W1 = rbinom(I, size = 1, prob=0.5) #binary
  
  #multinominal, prob={(1 − p)2, 2p(1 − p), p2}, where p = 40% 
  W = rmultinom(I,c(1,2,3),prob = c(0.36,0.48,0.16))
  for (obj in 1:I) {
    W2[obj]=which(W[,obj]==1)
  }
  W3 <- rnorm(I) #continuous
  dat <- data.frame(W1,W2,W3)
  colnames(dat) <- paste0("W",1:3)
  eta.surv = c(gamma%*%t(dat)) + c(alpha%*%t(eta.long[!duplicated(id),]))
  
  # simulate survival time
  scale = exp(-7)
  S = runif(I)
  Ti = rep(NA, I)
  alpha.beta = alpha%*%betat
  f = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c(alpha.beta)*t)
    }
    S[i] - exp(-stats::integrate(h, 0, tau)$value)
  }
  f = Vectorize(f)
  #curve(f,from = 0,to =100)
  
  for(i in 1:I){
    Ti[i] = uniroot(f, c(0, 100))$root
  }
  
  # simulate true survival probability
  pre.surtime = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c(alpha.beta)*t)
    }
    exp(-stats::integrate(h, 0, tau)$value)
  }
  
  OBS = obstime*12
  true.prob = matrix(NA, nrow=I, ncol=length(OBS[-1]))
  for(i in 1:I){
    ith = 0
    for(tau in OBS[-1]){
      ith = ith + 1
      true.prob[i, ith] = pre.surtime(tau)
    }
  }
  
  colnames(true.prob) = as.character(obstime[-1])
  
  Ti = Ti/12
  
  #--------------------------------
  # simulate censor time
  Left = Right = event = AR = AL = NA
  #######Interval censor
  for (i in 1:I) {
    A = cumsum(rexp(length(obstime),rate=1/MU)) #generate follow-up times
    A[1] = 0
    if(Ti[i]>max(A)) {
      AL[i]  <- max(A)
      AR[i] <- Inf}
    else {
      AL[i]  <- max(A[Ti[i]>A])
      AR[i] <- min(A[Ti[i]<=A])
    }
    
    if(is.infinite(AR[i])){
      Right[i] = Inf
      event[i] = 0
      if(AL[i]>max(obstime)){
        Left[i] = max(obstime)
      }else{
        Left[i] = max(obstime[AL[i]>obstime])
      }
    }else if(AL[i]>max(obstime)|AR[i]>max(obstime)){
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = Inf
      event[i] = 0
    }
    else {
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = min(obstime[AR[i]<obstime])
      event[i] = 1
    }
  }
  # prepare data
  visit = rep(1:J, I)
  obstime = rep(obstime, I) 
  erro = mvrnorm(N, c(0,0,0,0,0), diag(e.var))
  Y = matrix(0, nrow=N, ncol=5)
  for(i in 1:5)
    Y[,i] = eta.long[,i] + betat[i]*OBS + erro[,i]
  
  
  long.all = data.frame(id=id, visit=visit, time = rep(Ti, each=J), event = rep(event, each=J),
                        Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],Y4=Y[,4],Y5=Y[,5],obstime=obstime, X=X, ranef=ranef, dat[rep(1:I, each=J),], erro=I(erro))
  
  long = long.all
  
  # introduce missing complete at random
  if(miss){
    miss.index = sample(which(long$obstime>obstime[2]), 0.1*N)
    long = long[!c(1:N)%in%miss.index, ]
  }
  
  
  surv = data.frame(id = c(1:I),time=Ti, event=event,Left = Left, Right = Right, true.prob=I(true.prob), dat)
  
  # remove observations after event or censoring
  long = long[long$obstime<long$time, ]
  
  return(list(long=long, surv=surv, long.all=long.all))
}

sim_multi_nonlinear = function(I, J, MU, obstime = obstime, miss = FALSE, miss.rate = 0.1){
  # I : number of subjects
  # J : number of visits
  # obstime: observation times
  # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
  # miss.rate: missing rate.
  
  N = I*J
  # longitudinal submodel  
  beta0 = c(1,1,1,2,1)
  # beta1 = c(1.5,1,-1,1.5,1)
  # betat = c(0.5,0.5,1.2,-1, 1.5) #初始
  beta1 = c(0.5,1,-1,1.5,1)
  betat = c(1.5,1.5,1.2,-1.5, 1.5)
  ####generate subject-specific random effects from multi-normal distribution
  b.eta = runif(10,-0.5,0.5)
  b.sigma2 = c(0.5,1,1.5,2,2.5)
  e.var = c(1,1,1,1,1)
  bSigma = diag(b.sigma2)
  bSigma[1,2] = bSigma[2,1] = sqrt(b.sigma2[1]*b.sigma2[2])*b.eta[1]
  bSigma[1,3] = bSigma[3,1] = sqrt(b.sigma2[1]*b.sigma2[3])*b.eta[2]
  bSigma[1,4] = bSigma[4,1] = sqrt(b.sigma2[1]*b.sigma2[4])*b.eta[3]
  bSigma[1,5] = bSigma[5,1] = sqrt(b.sigma2[1]*b.sigma2[5])*b.eta[4]
  bSigma[2,3] = bSigma[3,2] = sqrt(b.sigma2[2]*b.sigma2[3])*b.eta[5]
  bSigma[2,4] = bSigma[4,2] = sqrt(b.sigma2[2]*b.sigma2[4])*b.eta[6]
  bSigma[2,5] = bSigma[5,2] = sqrt(b.sigma2[2]*b.sigma2[5])*b.eta[7]
  bSigma[3,4] = bSigma[4,3] = sqrt(b.sigma2[3]*b.sigma2[4])*b.eta[8]
  bSigma[3,5] = bSigma[5,3] = sqrt(b.sigma2[3]*b.sigma2[5])*b.eta[9]
  bSigma[4,5] = bSigma[5,4] = sqrt(b.sigma2[4]*b.sigma2[5])*b.eta[10]
  
  
  # sample covariate
  X = rep(rnorm(I, 6, 1), each=J)
  # sample random effect
  ranef = mvrnorm(I, c(0,0,0,0,0), bSigma)
  id = rep(1:I,each=J)
  ranef = ranef[id,]
  # construct longitudinal submodel
  eta.long = matrix(0, nrow=N, ncol=5)
  for(i in 1:5){
    eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
  }
  
  
  
  # survival submodel
  # alpha = c(0.5, -0.6, 0.6, 0.4, 0.2) # coefficients of time-varying
  # gamma = c(1.2, -0.8,-0.4) # coefficients of time-fixed
  alpha = c(0.5, -0.6, 0.6, 0.4, 0.2) # coefficients of time-varying
  gamma = c(1.2, -1.4,-0.4) # coefficients of time-fixed
  W1 = W2 = W3 = W = NA
  W1 = rbinom(I, size = 1, prob=0.5) #binary
  
  #multinominal, prob={(1 − p)2, 2p(1 − p), p2}, where p = 40% 
  W = rmultinom(I,c(1,2,3),prob = c(0.36,0.48,0.16))
  for (obj in 1:I) {
    W2[obj]=which(W[,obj]==1)
  }
  W3 <- rnorm(I) #continuous
  dat <- data.frame(W1,W2,W3)
  colnames(dat) <- paste0("W",1:3)
  eta.surv = c(gamma%*%t(dat)) + c(alpha%*%t(eta.long[!duplicated(id),]))
  
  # simulate survival time
  scale = exp(-7)
  S = runif(I)
  Ti = rep(NA, I)
  alpha.beta = alpha%*%betat
  
  c = c(1.2, 0.7, 0.4)
  knots = c(10, 20)
  f = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
                   c[3]*c(alpha.beta)* ifelse(t>=knots[2], t-knots[2], 0) )
    }
    S[i] - exp(-stats::integrate(h, 0, tau)$value)
  }
  f = Vectorize(f)
  curve(f,from = 0,to =100)
  
  
  for(i in 1:I){
    Ti[i] = uniroot(f, c(0,100))$root
  }
  sum(Ti/12<1)
  sum(Ti/12<2)
  sum(Ti/12<3)
  sum(Ti/12<4)
  
  # simulate true survival probability
  pre.surtime = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
                   c[3]*c(alpha.beta)* ifelse(t>=knots[2], t-knots[2], 0) ) 
      
    }
    exp(-stats::integrate(h, 0, tau)$value)
  }
  
  OBS = obstime*12
  true.prob = matrix(NA, nrow=I, ncol=length(OBS[-1]))
  for(i in 1:I){
    ith = 0
    for(tau in OBS[-1]){
      ith = ith + 1
      true.prob[i, ith] = pre.surtime(tau)
    }
  }
  
  colnames(true.prob) = as.character(obstime[-1])
  
  Ti = Ti/12
  
  #--------------------------------
  # simulate censor time
  Left = Right = event = AR = AL = NA
  #######Interval censor
  for (i in 1:I) {
    A = cumsum(rexp(length(obstime),rate=1/MU)) #generate follow-up times
    A[1] = 0
    #u = runif(length(obstime))
    # ind = ifelse(u<=p,1,0)
    # A = a[ind==1]
    # obtain intervals for dat
    # if(Ti[i]>max(A)) { #right censor, event=0
    #   Left[i]  <- max(A)
    #   Right[i] <- Inf
    #   event[i] = 0
    # } else if (Ti[i]<=min(A)) { #left censor, event=2
    #   Left[i] <- 0
    #   Right[i] <- min(A)
    #   event[i] = 1
    # } else {
    #   Left[i]  <- max(A[Ti[i]>A])
    #   Right[i] <- min(A[Ti[i]<=A])
    #   event[i] = 1
    # }
    
    if(Ti[i]>max(A)) {
      AL[i]  <- max(A)
      AR[i] <- Inf}
    else {
      AL[i]  <- max(A[Ti[i]>A])
      AR[i] <- min(A[Ti[i]<=A])
    }
    
    if(is.infinite(AR[i])){
      Right[i] = Inf
      event[i] = 0
      if(AL[i]>max(obstime)){
        Left[i] = max(obstime)
      }else{
        Left[i] = max(obstime[AL[i]>obstime])
      }
    }else if(AL[i]>max(obstime)|AR[i]>max(obstime)){
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = Inf
      event[i] = 0
    }
    else {
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = min(obstime[AR[i]<obstime])
      event[i] = 1
    }
  }
  
  
  
  # prepare data
  visit = rep(1:J, I)
  obstime = rep(obstime, I) 
  erro = mvrnorm(N, c(0,0,0,0,0), diag(e.var))
  Y = matrix(0, nrow=N, ncol=5)
  for(i in 1:5)
    Y[,i] = eta.long[,i] + betat[i]*(c[1]*ifelse(OBS<knots[1], OBS, knots[1]) + c[2]*ifelse(OBS<knots[2], pmax(0,(OBS-knots[1])), knots[2]-knots[1]) +
                                       c[3]* ifelse(OBS>=knots[2], OBS-knots[2], 0))  + erro[,i]
  
  
  long.all = data.frame(id=id, visit=visit, time = rep(Ti, each=J), event = rep(event, each=J),
                        Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],Y4=Y[,4],Y5=Y[,5],obstime=obstime, X=X, ranef=ranef,dat[rep(1:I, each=J),], erro=I(erro))
  
  long = long.all
  
  # introduce missing complete at random
  if(miss){
    miss.index = sample(which(long$obstime>obstime[2]), 0.1*N)
    long = long[!c(1:N)%in%miss.index, ]
  }
  
  surv = data.frame(id = c(1:I),time=Ti, event=event, Left = Left, Right = Right, true.prob=I(true.prob), dat)
  # remove observations after event or censoring
  long = long[long$obstime<long$time, ]
  
  return(list(long=long, surv=surv, long.all=long.all))
}

sim_ht_nonlinear = function(I, J, MU, obstime = obstime, miss = FALSE, miss.rate = 0.1){
  
  # I : number of subjects
  # J : number of visits
  # obstime: observation times
  # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
  # miss.rate: missing rate.
  
  N = I*J
  
  # longitudinal submodel  
  beta0 =c(1,0.5,1,2,1)
  beta1 = c(0.5,0.5,-4,0.5,1)
  betat = c(1.5,1.5,2,-3, 1.5) #初始
  # beta1 = c(1.5,1,-1,1.5,1)
  # betat = c(1.5,1.5,1.2,-1.5, 1.5)
  ####generate subject-specific random effects from multi-normal distribution
  b.eta = runif(10,-0.5,0.5)
  b.sigma2 = b.sigma2 = c(0.5,1,1.5,2,2.5)
  e.var = c(1,1,1,1,1)
  bSigma = diag(b.sigma2)
  bSigma[1,2] = bSigma[2,1] = sqrt(b.sigma2[1]*b.sigma2[2])*b.eta[1]
  bSigma[1,3] = bSigma[3,1] = sqrt(b.sigma2[1]*b.sigma2[3])*b.eta[2]
  bSigma[1,4] = bSigma[4,1] = sqrt(b.sigma2[1]*b.sigma2[4])*b.eta[3]
  bSigma[1,5] = bSigma[5,1] = sqrt(b.sigma2[1]*b.sigma2[5])*b.eta[4]
  bSigma[2,3] = bSigma[3,2] = sqrt(b.sigma2[2]*b.sigma2[3])*b.eta[5]
  bSigma[2,4] = bSigma[4,2] = sqrt(b.sigma2[2]*b.sigma2[4])*b.eta[6]
  bSigma[2,5] = bSigma[5,2] = sqrt(b.sigma2[2]*b.sigma2[5])*b.eta[7]
  bSigma[3,4] = bSigma[4,3] = sqrt(b.sigma2[3]*b.sigma2[4])*b.eta[8]
  bSigma[3,5] = bSigma[5,3] = sqrt(b.sigma2[3]*b.sigma2[5])*b.eta[9]
  bSigma[4,5] = bSigma[5,4] = sqrt(b.sigma2[4]*b.sigma2[5])*b.eta[10]
  
  
  # sample covariate
  X = rep(rnorm(I, 3, 1), each=J)
  # sample random effect
  ranef = mvrnorm(I, c(0,0,0,0,0), bSigma)
  id = rep(1:I,each=J)
  ranef = ranef[id,]
  # construct longitudinal submodel
  eta.long = matrix(0, nrow=N, ncol=5)
  for(i in 1:5){
    eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
  }
  
  
  
  #survival submodel
  alpha = c(0.35, -0.6, 0.6, 0.2, 0.1) # coefficients of time-varying
  gamma = c(0.2, -0.1,-0.2,4) # coefficients of time-fixed
  
  
  W1 = W2 = W3 = W4 = W5 = W = dat = NA
  W1 = rbinom(I, size = 1, prob=0.5) #binary
  
  #multinominal, prob={(1 − p)2, 2p(1 − p), p2}, where p = 40% 
  W = rmultinom(I,c(1,2,3),prob = c(0.36,0.48,0.16))
  for (obj in 1:I) {
    W2[obj]=which(W[,obj]==1)
  }
  W3 <- rnorm(I) #continuous
  W4 <- c(W1*W3) #interaction
  #W5 <- c(W1*W1) 
  dat <- data.frame(W1,W2,W3,W4)
  colnames(dat) <- paste0("W",1:4)
  eta.surv = c(gamma%*%t(dat)) + c(alpha%*%t(eta.long[!duplicated(id),]))
  
  # simulate survival time
  scale = exp(-7)
  S = runif(I)
  Ti = rep(NA, I)
  alpha.beta = alpha%*%betat
  f = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c(alpha.beta)*t)
    }
    S[i] - exp(-stats::integrate(h, 0, tau)$value)
  }
  f = Vectorize(f)
  curve(f,from = 0,to =100)
  
  for(i in 1:I){
    Ti[i] = uniroot(f, c(0, 100))$root
  }  
  
  sum(Ti/12<1)
  sum(Ti/12<2)
  sum(Ti/12<3)
  sum(Ti/12<4)
  
  # simulate true survival probability
  pre.surtime = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c(alpha.beta)*t)
    }
    exp(-stats::integrate(h, 0, tau)$value)
  }
  
  OBS = obstime*12
  true.prob = matrix(NA, nrow=I, ncol=length(OBS[-1]))
  for(i in 1:I){
    ith = 0
    for(tau in OBS[-1]){
      ith = ith + 1
      true.prob[i, ith] = pre.surtime(tau)
    }
  }
  
  colnames(true.prob) = as.character(obstime[-1])
  
  Ti = Ti/12
  
  #--------------------------------
  # simulate censor time
  Left = Right = event = AR = AL = NA
  #######Interval censor
  for (i in 1:I) {
    A = cumsum(rexp(length(obstime),rate=1/MU)) #generate follow-up times
    A[1] = 0
    if(Ti[i]>max(A)) {
      AL[i]  <- max(A)
      AR[i] <- Inf}
    else {
      AL[i]  <- max(A[Ti[i]>A])
      AR[i] <- min(A[Ti[i]<=A])
    }
    
    if(is.infinite(AR[i])){
      Right[i] = Inf
      event[i] = 0
      if(AL[i]>max(obstime)){
        Left[i] = max(obstime)
      }else{
        Left[i] = max(obstime[AL[i]>obstime])
      }
    }else if(AL[i]>max(obstime)|AR[i]>max(obstime)){
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = Inf
      event[i] = 0
    }
    else {
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = min(obstime[AR[i]<obstime])
      event[i] = 1
    }
  }
  # prepare data
  visit = rep(1:J, I)
  obstime = rep(obstime, I) 
  erro = mvrnorm(N, c(0,0,0,0,0), diag(e.var))
  Y = matrix(0, nrow=N, ncol=5)
  for(i in 1:5)
    Y[,i] = eta.long[,i] + betat[i]*OBS + erro[,i]
  
  
  long.all = data.frame(id=id, visit=visit, time = rep(Ti, each=J), event = rep(event, each=J),
                        Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],Y4=Y[,4],Y5=Y[,5],obstime=obstime, X=X, ranef=ranef, dat[rep(1:I, each=J),], erro=I(erro))
  
  long = long.all
  
  # introduce missing complete at random
  if(miss){
    miss.index = sample(which(long$obstime>obstime[2]), 0.1*N)
    long = long[!c(1:N)%in%miss.index, ]
  }
  
  
  surv = data.frame(id = c(1:I),time=Ti, event=event,Left = Left, Right = Right, true.prob=I(true.prob), dat)
  
  # remove observations after event or censoring
  long = long[long$obstime<long$time, ]
  
  return(list(long=long, surv=surv, long.all=long.all))
}

sim_both_nonlinear = function(I, J, MU, obstime = obstime, miss = FALSE, miss.rate = 0.1){
  # I : number of subjects
  # J : number of visits
  # obstime: observation times
  # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
  # miss.rate: missing rate.
  
  N = I*J
  # longitudinal submodel  
  beta0 =c(1,0.5,1,2,1)
  beta1 = c(0.5,0.5,-4,0.5,1)
  betat = c(1.5,1.5,2,-3, 1.5) #初始
  ####generate subject-specific random effects from multi-normal distribution
  b.eta = runif(10,-0.5,0.5)
  b.sigma2 = c(0.5,1,1.5,2,2.5)
  e.var = c(1,1,1,1,1)
  bSigma = diag(b.sigma2)
  bSigma[1,2] = bSigma[2,1] = sqrt(b.sigma2[1]*b.sigma2[2])*b.eta[1]
  bSigma[1,3] = bSigma[3,1] = sqrt(b.sigma2[1]*b.sigma2[3])*b.eta[2]
  bSigma[1,4] = bSigma[4,1] = sqrt(b.sigma2[1]*b.sigma2[4])*b.eta[3]
  bSigma[1,5] = bSigma[5,1] = sqrt(b.sigma2[1]*b.sigma2[5])*b.eta[4]
  bSigma[2,3] = bSigma[3,2] = sqrt(b.sigma2[2]*b.sigma2[3])*b.eta[5]
  bSigma[2,4] = bSigma[4,2] = sqrt(b.sigma2[2]*b.sigma2[4])*b.eta[6]
  bSigma[2,5] = bSigma[5,2] = sqrt(b.sigma2[2]*b.sigma2[5])*b.eta[7]
  bSigma[3,4] = bSigma[4,3] = sqrt(b.sigma2[3]*b.sigma2[4])*b.eta[8]
  bSigma[3,5] = bSigma[5,3] = sqrt(b.sigma2[3]*b.sigma2[5])*b.eta[9]
  bSigma[4,5] = bSigma[5,4] = sqrt(b.sigma2[4]*b.sigma2[5])*b.eta[10]
  
  
  # sample covariate
  X = rep(rnorm(I, 6, 1), each=J)
  # sample random effect
  ranef = mvrnorm(I, c(0,0,0,0,0), bSigma)
  id = rep(1:I,each=J)
  ranef = ranef[id,]
  # construct longitudinal submodel
  eta.long = matrix(0, nrow=N, ncol=5)
  for(i in 1:5){
    eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
  }
  
  
  
  # survival submodel
  alpha = c(0.8, -0.6, 0.6, 0.2, 0.1) # coefficients of time-varying
  gamma = c(0.2, -0.1,-0.2,4) # coefficients of time-fixed
  
  
  W1 = W2 = W3 = W4 = W5 = W = dat = NA
  W1 = rbinom(I, size = 1, prob=0.5) #binary
  
  #multinominal, prob={(1 − p)2, 2p(1 − p), p2}, where p = 40% 
  W = rmultinom(I,c(1,2,3),prob = c(0.36,0.48,0.16))
  for (obj in 1:I) {
    W2[obj]=which(W[,obj]==1)
  }
  W3 <- rnorm(I) #continuous
  W4 <- c(W1*W3) #interaction
  #W5 <- c(W1*W1) 
  dat <- data.frame(W1,W2,W3,W4)
  colnames(dat) <- paste0("W",1:4)
  eta.surv = c(gamma%*%t(dat)) + c(alpha%*%t(eta.long[!duplicated(id),]))
  
  # simulate survival time
  scale = exp(-7)
  S = runif(I)
  Ti = rep(NA, I)
  alpha.beta = alpha%*%betat
  # c = c(1, 0.7, 0.5)
  # knots = c(6, 15,25)
  # f = function(tau){
  #   h = function(t) {
  #     scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
  #                  c[3]*c(alpha.beta)* ifelse(t>= knots[2], t-knots[2], 0) )
  #   }
  #   S[i] - exp(-stats::integrate(h, 0, tau)$value)
  # }
  c = c(1.2, 0.7, 0.5)
  knots = c(8, 20)
  f = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
                   c[3]*c(alpha.beta)* ifelse(t>=knots[2], t-knots[2], 0))
    }
    S[i] - exp(-stats::integrate(h, 0, tau)$value)
  }
  f = Vectorize(f)
  curve(f,from = 0,to =100)
  
  for(i in 1:I){
    Ti[i] = uniroot(f, c(0,100))$root
  }
  
  # sum(Ti/12<0.5)
  # sum(Ti/12<1)
  # sum(Ti/12<1.5)
  # sum(Ti/12<2)
  # sum(Ti/12<2.5)
  # simulate true survival probability
  pre.surtime = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
                   c[3]*c(alpha.beta)* ifelse(t>=knots[2], t-knots[2], 0)) 
      
    }
    exp(-stats::integrate(h, 0, tau)$value)
  }
  
  OBS = obstime*12
  true.prob = matrix(NA, nrow=I, ncol=length(OBS[-1]))
  for(i in 1:I){
    ith = 0
    for(tau in OBS[-1]){
      ith = ith + 1
      true.prob[i, ith] = pre.surtime(tau)
    }
  }
  
  colnames(true.prob) = as.character(obstime[-1])
  
  Ti = Ti/12
  
  #--------------------------------
  # simulate censor time
  Left = Right = event = AR = AL = NA
  #######Interval censor
  for (i in 1:I) {
    A = cumsum(rexp(length(obstime),rate=1/MU)) #generate follow-up times
    A[1] = 0
    #u = runif(length(obstime))
    # ind = ifelse(u<=p,1,0)
    # A = a[ind==1]
    # obtain intervals for dat
    # if(Ti[i]>max(A)) { #right censor, event=0
    #   Left[i]  <- max(A)
    #   Right[i] <- Inf
    #   event[i] = 0
    # } else if (Ti[i]<=min(A)) { #left censor, event=2
    #   Left[i] <- 0
    #   Right[i] <- min(A)
    #   event[i] = 1
    # } else {
    #   Left[i]  <- max(A[Ti[i]>A])
    #   Right[i] <- min(A[Ti[i]<=A])
    #   event[i] = 1
    # }
    
    if(Ti[i]>max(A)) {
      AL[i]  <- max(A)
      AR[i] <- Inf}
    else {
      AL[i]  <- max(A[Ti[i]>A])
      AR[i] <- min(A[Ti[i]<=A])
    }
    
    if(is.infinite(AR[i])){
      Right[i] = Inf
      event[i] = 0
      if(AL[i]>max(obstime)){
        Left[i] = max(obstime)
      }else{
        Left[i] = max(obstime[AL[i]>obstime])
      }
    }else if(AL[i]>max(obstime)|AR[i]>max(obstime)){
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = Inf
      event[i] = 0
    }
    else {
      Left[i] = max(obstime[AL[i]>=obstime])
      Right[i] = min(obstime[AR[i]<obstime])
      event[i] = 1
    }
  }
  
  
  
  # prepare data
  visit = rep(1:J, I)
  obstime = rep(obstime, I) 
  erro = mvrnorm(N, c(0,0,0,0,0), diag(e.var))
  Y = matrix(0, nrow=N, ncol=5)
  for(i in 1:5)
    Y[,i] = eta.long[,i] + betat[i]*(c[1]*ifelse(OBS<knots[1], OBS, knots[1]) + c[2]*ifelse(OBS<knots[2], pmax(0,(OBS-knots[1])), knots[2]-knots[1]) +
                                       c[3]* ifelse(OBS>=knots[2], OBS-knots[2], 0))  + erro[,i]
  
  
  long.all = data.frame(id=id, visit=visit, time = rep(Ti, each=J), event = rep(event, each=J),
                        Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],Y4=Y[,4],Y5=Y[,5],obstime=obstime, X=X, ranef=ranef,dat[rep(1:I, each=J),], erro=I(erro))
  
  long = long.all
  
  # introduce missing complete at random
  if(miss){
    miss.index = sample(which(long$obstime>obstime[2]), 0.1*N)
    long = long[!c(1:N)%in%miss.index, ]
  }
  
  surv = data.frame(id = c(1:I),time=Ti, event=event, Left = Left, Right = Right, true.prob=I(true.prob), dat)
  # remove observations after event or censoring
  long = long[long$obstime<long$time, ]
  
  return(list(long=long, surv=surv, long.all=long.all))
}
