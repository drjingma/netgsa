glmnet.soft <-
function(x, y, lambda){
  ## Now takes multiple values of lambda
  betaMat <- matrix(0, nrow = ncol(x), ncol = length(lambda))
  
  ##If the structure is certain (weight can be set to zero), then simply return 1's
  betaMat[,which(lambda == 0)] <- 1 
  
  non0_lambdas <- lambda[lambda != 0]
  any_non0 <- length(non0_lambdas) > 0
  if(any_non0){
    if(ncol(x) > 1){		## use glmnet
      fit = glmnet::glmnet(x, y, family="gaussian", alpha=1, lambda=non0_lambdas)
      beta = as.matrix(fit$beta)
    }else{				    ## use soft thresholding
      beta = lm(y~x)$coef[2]
      beta = matrix(sign(beta) * pmax((abs(beta) - non0_lambdas/2),0), nrow = 1) ##pmax so we can have multiple lambdas
    } 
    betaMat[,which(lambda == non0_lambdas)] <- beta
  }
  
  return(betaMat)
}
