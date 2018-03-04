glmnet.soft <-
function(x, y, lambda){
  if (lambda==0){
    ##If the structure is certain (weight can be set to zero), then simply return 1's
    p = ncol(matrix(x, nrow=length(y)))
    beta = rep(1, p)
  } else {
    if(ncol(x) > 1){		## use glmnet
      fit = glmnet::glmnet(x, y, family="gaussian", alpha=1, lambda=lambda)
      beta = as.matrix(fit$beta)
    }else{				    ## use soft thresholding
      beta = lm(y~x)$coef[2]
      beta = sign(beta) * max((abs(beta) - lambda/2),0)
    }
  }
  return(beta)
}
