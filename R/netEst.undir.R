netEst.undir <-
  function(x, zero = NULL, one = NULL, lambda, rho=NULL, penalize_diag = TRUE, weight = NULL, eta=0, verbose = FALSE, eps = 1e-08) {
    p <- nrow(x)
    n <- ncol(x)
    
    Adj = lapply(1:length(lambda), function(x) matrix(0, p, p))
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
    
    if (any(abs(lambda) < eps)){
      stop("The penalty parameter lambda needs to be greater than zero!")
    }
    
    if (n<10){
      warning("The sample size is too small! Network estimate may be unreliable!")
    }
    
    if (is.null(weight)){
      weight = 0
    }
    
    if (is.null(rho)) {
      rhoM = matrix(0.1*sqrt(log(p)/n), p, p)
    } else if (is.matrix(rho)){
      if(length(rho)!=p*p) 
        stop("The input matrix for \"rho\" must be of size ",p," by ",p)
      rhoM = rho             
    } else {
      rhoM = matrix(rho,p,p) 
    }
    
    X = t(x)
    X = scale(X, center = TRUE, scale = TRUE)
    
    # check if network information is complete
    if ( identical(one+zero, matrix(1,p,p)-diag(p)) && (weight == 0)){
      cat("Network information is complete! \n")
      Adj = lapply(Adj, function(nothing) return(one))
    } else {
      if (!is.null(weight) ){
        if (weight < -1e-16){
          stop("Negative weight parameter detected! Please double check!")
        }
      }
      
      # estimate the network topology
      for (i in 1:p) {
        Y = matrix(X[, i], ncol = 1)
        Xmat = X[, -i, drop = FALSE]
        
        ## Get the zero and one indices. 
        infoInd = one[i, -i] - zero[i, -i]
        beta = matrix(0, p - 1, length(lambda)) ## For multiple lambda sequence
        
        if (sum(infoInd == 0) == 0) {
          if (verbose) {
            cat("Complete information known! \n ")
          }
          if (sum(infoInd == 1)>0){
            Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
            beta[(infoInd == 1), ] = glmnet.soft(Xmat1, Y, lambda = lambda*weight)
          }
          if (sum(infoInd ==-1)>0){
            beta[(infoInd ==-1), ] = 0              
          }
        } else {
          if (verbose) {
            cat("Incomplete information known! \n ")
          }
          if (sum(infoInd == -1) == 0) {
            if (sum(infoInd == 1) == 0) {
              if (verbose) {
                cat("Incomplete information: no 0's and no 1's! \n ")
              }
              beta = glmnet.soft(Xmat, Y, lambda = lambda)
            } else {
              if (verbose) {
                cat("Incomplete information: no 0's, but with 1's! \n ")
              }
              if (sum(infoInd==1)>=n) {#if there are fewer samples than the number of 1's, simply adopt the one's and pass
                beta[which(infoInd == 1), ] = 1 
                beta[which(infoInd == 0), ] = 0 
              } else {
                Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
                Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
                beta[(infoInd == 1), ] = glmnet.soft(Xmat1, Y, lambda = lambda*weight)
                tmp = as.matrix(beta[(infoInd == 1), , drop = FALSE])
                tmp_col_fit = do.call(cbind,
                                      lapply(1:length(lambda), function(i){
                                        res.glm = Y - Xmat1 %*Cpp% tmp[,i]
                                        return(glmnet.soft(Xmat2, res.glm, lambda = lambda[i]))
                                      }))
                beta[(infoInd == 0), ] = tmp_col_fit
              }
            }
          } else {
            if (sum(infoInd == 1) == 0) {
              if (verbose) {
                cat("Incomplete information: with 0's and no 1's! \n ")
              }
              beta[(infoInd == -1), ] = 0
              Xnew = Xmat[, (infoInd != -1), drop = FALSE]
              beta[(infoInd != -1), ] = glmnet.soft(Xnew, Y, lambda = lambda)
            } else {
              if (verbose) {
                cat("Incomplete information: with both 0's and 1's! \n ")
              }
              beta[(infoInd == -1), ] = 0 #known non-edges 
              Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
              Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
              beta[(infoInd == 1), ] = glmnet.soft(Xmat1, Y, lambda = lambda*weight) 
              tmp = as.matrix(beta[(infoInd == 1), , drop = FALSE])
              tmp_col_fit = do.call(cbind,
                                    lapply(1:length(lambda), function(i){
                                          res.glm = Y - Xmat1 %*Cpp% tmp[,i]
                                          return(glmnet.soft(Xmat2, res.glm, lambda = lambda[i]))
                                      }))
              beta[(infoInd == 0), ] = tmp_col_fit
            }
          }
        }
        # Fill Adj[[1]] row with beta[,1], and Adj[[2]] with beta[,2], etc
        Adj <- Map(function(Amat, lam_num)  {Amat[i, -i] = as.vector(beta[,lam_num])
                                             return(Amat)}, Adj, 1:length(lambda))
      }
      
      ## symmetrization 
      Adj = lapply(Adj, function(Amat){
        temp = (Amat + t(Amat))/2
        temp = (abs(temp) > eps)
        diag(temp) = 0
        return(temp)
        })
    }

    empcov = cov(X)
    if (kappa(empcov) > 1e+3){
      empcov = empcov + eta * diag(p)
    }
    diag_penalty <- ifelse(penalize_diag, 0.01*sqrt(log(p)/n), 0)
    BIG = 10e9
    partialCor_l <- lapply(Adj, function(Amat){
      rhoM[which(Amat==0)] = BIG;
      diag(rhoM) = diag_penalty
      
      #No need to handle singleton clusters. They shouldn't exist here. Are handled in netEstClusts function
      obj <- glassoFast(empcov, rho = rhoM, thr = 1e-4)
      siginv <- chol2inv(cholCpp(obj$w))

      partialCor <- Ip - cov2cor(siginv)
      partialCor[abs(partialCor) < 1e-08] <- 0
      rownames(partialCor) <- rownames(x); colnames(partialCor) <- rownames(x);
      rownames(siginv) <- rownames(x); colnames(siginv) <- rownames(x);
      out <- list(partialCor=partialCor,siginv=siginv,lambda=lambda)
      #out <- calcPartialCor(Amat, rhoM, empcov, rownames(x))
      return(out)
    })

    partialCor <- lapply(partialCor_l, "[[", "partialCor")
    siginv     <- lapply(partialCor_l, "[[", "siginv")


    return(list(Adj=partialCor,invcov=siginv,lambda=lambda))
  }
