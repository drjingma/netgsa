netEst.undir <-
function(X, zero = NULL, one = NULL, lambda, rho = NULL, weight= NULL, eta=0, verbose = FALSE, eps = 1e-08) {
  n = dim(X)[1]
  p = dim(X)[2]
  Adj = matrix(0, p, p)
  Ip = diag(rep(1, p))
  
  if (is.null(zero)) {
    zero = matrix(0, p, p)
  }
  
  if (is.null(one)) {
    one = matrix(0, p, p)
  }
  
  if (abs(lambda) < eps){
    stop("The penalty parameter lambda needs to be greater than zero!")
  }
  
  if (is.null(rho)) {
    rho = 0.1*sqrt(log(p)/n)
  }
  
  if (!is.null(weight) ){
    if (weight < -1e-16){
      stop("Negative weight parameter detected! Please double check!")
    }
  }
  
  if (is.null(weight)){
    weight = 0
  }
  
  ## To get the empirical covariance matrix 
  X = scale(X, center = TRUE, scale = TRUE)
  for (i in 1:p) {
    # print(i)
    Y = matrix(X[, i], ncol = 1)
    Xmat = X[, -i]
    
    ## Get the zero and one indices. 
    infoInd = one[i, -i] - zero[i, -i]
    beta = matrix(0, p - 1, 1)
    
    if (sum((infoInd == 0)) == 0) {
      if (verbose) {
        cat("Complete information known!  ")
      }
      Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
      beta[(infoInd == 1), ] = glmnet.soft(Xmat1, Y, lambda = lambda*weight)
      beta[(infoInd == 0), ] = 0              
    } else {
      if (verbose) {
        cat("Incomplete information known!  ")
      }
      
      if (sum((infoInd == -1)) == 0) {
        if (sum((infoInd == 1)) == 0) {
          if (verbose) {
            cat("Incomplete information: no 0's and no 1's!  ")
          }
          beta = glmnet.soft(Xmat, Y, lambda = lambda)
        } else {
          if (verbose) {
            cat("Incomplete information: no 0's, but with 1's!  ")
          }
          if (sum(infoInd==1)>=n) {#if there are fewer samples than the number of 1's, simply adopt the one's and pass
            beta[which(infoInd == 1), ] = 1 
            beta[which(infoInd == 0), ] = 0 
          } else {
            Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
            Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
            beta[(infoInd == 1), ] = glmnet.soft(Xmat1, Y, lambda = lambda*weight)
            tmp = as.matrix(beta[(infoInd == 1), ])
            res.glm = Y - Xmat1 %*% tmp
            beta[(infoInd == 0), ] = glmnet.soft(Xmat2, res.glm, lambda = lambda)
          }
        }
      } else {
        if (sum((infoInd == 1)) == 0) {
          if (verbose) {
            cat("Incomplete information: with 0's and no 1's!  ")
          }
          beta[(infoInd == -1), ] = 0
          Xnew = Xmat[, (infoInd != -1)]
          beta[(infoInd != -1), ] = glmnet.soft(Xnew, Y, lambda = lambda)
        } else {
          if (verbose) {
            cat("Incomplete information: with both 0's and 1's!  ")
          }
          beta[(infoInd == -1), ] = 0 #known non-edges 
          Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
          Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
          beta[(infoInd == 1), ] = glmnet.soft(Xmat1, Y, lambda = lambda*weight) 
          tmp = as.matrix(beta[(infoInd == 1), ])
          res.glm = Y - Xmat1 %*% tmp
          beta[(infoInd == 0), ] = glmnet.soft(Xmat2, res.glm, lambda = lambda)
        }
      }
    }
    Adj[i, -i] = as.vector(beta);
  }
  
  ## Symmetrization 
  Adj = (Adj + t(Adj))/2
  Adj = (abs(Adj) > eps)
  diag(Adj) = 0; 

  ## Estimate the partial correlation matrix based on graphical lasso 
  info = zeroInd(Adj, 1)$zeroArr
  empcov <- cov(X) #empirical cov
  if (kappa(empcov) > 1e+3){
    empcov = empcov + eta * diag(p)
  }
  obj <- glasso(empcov, rho = rho, zero = info, penalize.diagonal = FALSE)
  siginv = chol2inv(chol(obj$w))
  
  partialCor = Ip - cov2cor(siginv)
  partialCor[abs(partialCor) < 1e-08] <- 0
  
  return(list(Adj = partialCor, invcov=siginv,lambda=lambda))
}
