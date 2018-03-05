netEst.dir <-
function(x, zero=NULL, one=NULL, lambda = NULL, verbose=FALSE, eps=1e-08) {
  p <- nrow(x)
  n <- ncol(x)
  
  Adj = matrix(0, p, p)
  Ip = diag(rep(1, p))
  
  if (is.null(zero)) {
    zero = matrix(0, p, p)
    rownames(zero) = rownames(x)
  } else {
    # check if variables match
    if (!identical(rownames(zero),rownames(x))){
      stop('Variables in the nonedge matrix and data do not match!')
    }
  }
  
  if (is.null(one)) {
    one = matrix(0, p, p)
    rownames(one) = rownames(x)
  }  else {
    # check if variables match
    if (!identical(rownames(one),rownames(x))){
      stop('Variables in the edge matrix and data do not match!')
    }
  }
  
  if (sum(one*zero) > 0){
    stop("Information on 0's and 1's overlaps!")
  }  
  if (n<10){
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  if (is.null(lambda)){
    alpha = 0.25
    lambda = 2 * qnorm(1-alpha/(2*p* seq(1,p-1) )) /sqrt(n)
  } else if (length(lambda)==1){
    lambda = c(0, rep(lambda, p-1))
  } else if (length(lambda) == p-1){
    lambda = c(0, lambda)
  } else {
    stop("The penalty parameter is incorrectly specified!")
  }
  
  ## To get the empirical correlation matrix 
  X <- t(x)
  X = scale(X, center=TRUE, scale=TRUE)
  for (i in 2:p) {
    Y = matrix(X[, i], ncol=1)
    Xmat = matrix(X[, 1:(i-1)], ncol=i-1)
    
    ## Get the zero and one indices 
    infoInd = one[i, 1:(i-1)] - zero[i, 1:(i-1)]
    beta = matrix(0, i-1, 1)
    
    if (sum((infoInd == 0)) == 0) {
      if (verbose) {
        cat("Complete information known!  ")
      }
      beta[which(infoInd == 1), ] = 1
      beta[which(infoInd == -1), ] = 0
    } else {
      if (verbose) {
        cat("Incomplete information known!  ")
      }
      
      if (sum((infoInd == -1)) == 0) {
        if (sum((infoInd == 1)) == 0) {
          if (verbose) {
            cat("Incomplete information: no 0's and no 1's!  ")
          }
          beta = glmnet.soft(x=Xmat, y=Y, lambda=lambda[i])	#glmnet if ncol(X) > 1, o.w. soft thresholding
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
            fit.lm = lm(Y ~ 0 + Xmat1)
            beta[(infoInd == 1), ] = as.numeric(coef(fit.lm))
            res.lm = as.numeric(residuals(fit.lm)) ##get the residuals           
            beta[(infoInd == 0), ] = glmnet.soft(x=Xmat2, y=res.lm, lambda=lambda[i])
          }
        }
      } else {
        if (sum((infoInd == 1)) == 0) {
          if (verbose) {
            cat("Incomplete information: with 0's and no 1's!  ")
          }
          beta[(infoInd == -1), ] = 0
          Xnew = Xmat[, (infoInd != -1)]
          beta[(infoInd != -1), ] = glmnet.soft(x=Xnew, y=Y, lambda=lambda[i])
        } else {
          if (verbose) {
            cat("Incomplete information: with both 0's and 1's!  ")
          }
          beta[(infoInd == -1), ] = 0 #known non-edges 
          Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
          Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
          fit.lm = lm(Y ~ 0 + Xmat1)
          beta[(infoInd == 1), ] = as.numeric(coef(fit.lm))
          res.lm = as.numeric(residuals(fit.lm)) ##get the residuals           
          beta[(infoInd == 0), ] = glmnet.soft(x=Xmat2, y=res.lm, lambda=lambda[i])
        }
      }
    }
    beta[(abs(beta) < eps)] = 0
    Adj[i, 1:(i-1)] = as.vector(beta); 
  }
  
  infmat = solve(Ip - Adj)
  infmat[abs(infmat) < eps] <- 0
  rownames(Adj) = rownames(x);
  colnames(Adj) = rownames(x);
  
  return(list(Adj=Adj, infmat=infmat, lambda=lambda))
}
