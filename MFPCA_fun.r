################### MFPCA Functions #######################

# uPACE       # univariate FPCA via PACE
# mfpca.score # MFPC scores
# mfpca.pred  # MFPC longitudinal prediction
# cond.prob   # Risk conditional probability
# minmax      # scale function


# univariate FPCA via principal analysis by conditional estimation(PACE)
uPACE = function(testData, domain, predData=NULL, nbasis = 10, pve = 0.95, npc = NULL){
  
  tmp = funData(domain, testData)
  if(is.null(predData)){
    tmp2 = NULL
  }else{
    tmp2 = funData(domain, predData)
  }
  
  res = PACE(tmp, tmp2, pve=pve, npc= npc, nbasis=nbasis)
  return(res)
} 

# multivariate FPCA based on results from uPACE
# mFPCA = function(Xi, phi, p , L){
#   
#   # eigenanalysis on matrix M
#   M = t(Xi)%*%Xi/(I-1)
#   eigen.M = eigen(M)
#   values = eigen.M$values
#   pve = cumsum(values)/sum(values)
#   Cms = eigen.M$vectors
#   index = unlist(lapply(1:length(L), function(x) rep(x, L[x])))
#   
#   # MFPCA score
#   rho = mfpca.score(Xi, Cms)
#   
#   # MFPCA eigenfunction
#   psis = NULL
#   for(j in 1:p){
#     psi = NULL
#     for(m in 1:dim(Cms)[2]){
#       psi = cbind(psi, phi[[j]]%*%Cms[which(index==j),m])
#     }
#     psis[[j]] = psi
#   }
#   
#   out = list(eigenvalue = values, Cms = Cms, pve = pve, index=index, rho = rho, psis=psis)
#   
#   return(out)
# }
# multivariate FPCA based on results from uPACE
mFPCA = function(Xi, phi, p , L, nsurv,predXi=NULL){
  
  # eigenanalysis on matrix M
  M = t(Xi)%*%Xi/(nsurv-1)
  eigen.M = eigen(M)
  values = eigen.M$values
  pve = cumsum(values)/sum(values)
  Cms = eigen.M$vectors
  index = unlist(lapply(1:length(L), function(x) rep(x, L[x])))
  
  # MFPCA score
  if(is.null(predXi)){
    predXi = Xi
  }
  rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
  for(i in 1:nrow(predXi)){
    for(m in 1:dim(Cms)[2]){
      rho[i,m] = predXi[i,]%*%Cms[,m]
    }
  }
  
  # MFPCA eigenfunction
  psis = NULL
  for(j in 1:p){
    psi = NULL
    for(m in 1:dim(Cms)[2]){
      psi = cbind(psi, phi[[j]]%*%Cms[which(index==j),m])
    }
    psis[[j]] = psi
  }
  
  out = list(eigenvalue = values, Cms = Cms, pve = pve, index=index, rho = rho, psis=psis)
  
  return(out)
}
# mfpc score calculation
mfpca.score = function(predXi, Cms){
  rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
  for(i in 1:nrow(predXi)){
    for(m in 1:dim(Cms)[2]){
      rho[i,m] = predXi[i,]%*%Cms[,m]
    }
  }
  return(rho)
}


# mfpc trajectories prediction
mfpca.pred = function(score, meanf, psi, n.rho=NULL){
  p = length(psi)
  n = nrow(score)
  
  if(is.null(n.rho)){
    n.rho = ncol(score)
  }
  
  pred = array(NA, c(n, length(meanf[[1]]), p))
  for(m in 1:p){
    pred[,,m] = matrix(meanf[[m]], nrow=n, ncol =length(meanf[[m]]), byrow = T ) + score[,1:n.rho]%*%t(psi[[m]][, 1:n.rho])
  }
  
  out = pred
  return(out)
}


#risk predict
cond.prob = function(model, newdata, Tstart, Tpred){
  risk.Tstart = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tstart)$surv)
  risk.Tpred = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tpred)$surv)
  return(risk.Tpred/risk.Tstart)
}

minmax <- function(x){
  newx <- (x - min(x[!is.na(x)]))/(max(x[!is.na(x)]) - min(x[!is.na(x)]))
  return(newx)
}