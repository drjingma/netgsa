bic.netEst.undir <-
function(x, zero = NULL, one = NULL, lambda, rho = NULL,weight = NULL, eta = 0, verbose = FALSE, eps = 1e-08) {
  p = nrow(x)
  n = ncol(x)
  Adj = array(0, c(p, p, length(lambda)))
  
  if (is.null(zero)) {
    zero = matrix(0, p, p)
    rownames(zero) = rownames(x)
  }
  
  if (is.null(one)) {
    one = matrix(0, p, p)
    rownames(one) = rownames(x)
  }
  
  if (min(lambda)==0){
    stop("The penalty parameter lambda needs to be greater than zero!")
  }
  
  if (!is.null(weight) ){
    if (weight < -1e-16){
      stop("Negative weight parameter detected! Please double check!")
    }
  }
  
  if (is.null(weight)){
    weight = 0
  }  

  bic.score = matrix(0, length(lambda), length(weight))
  df = matrix(0, length(lambda), length(weight))
  
  empcov = cov(t(x))
  
  if (kappa(empcov) > 1e+3) {
    empcov = empcov + eta * diag(p)
  }
  
  for (loop.lambda in 1:length(lambda)) {
    for (loop.w in 1:length(weight)){
      siginv = netEst.undir(x, zero, one, lambda[loop.lambda], rho=rho, weight=weight[loop.w], eta=eta)$invcov
      no.edge = sum(abs(siginv) > eps) - p
      bic.score[loop.lambda, loop.w] = matTr2(empcov %*% siginv) - determinant(siginv, logarithm = T)$modulus + log(n) * no.edge/(2 * n)
      df[loop.lambda, loop.w] = no.edge
    }
  }

  out <- list(lambda=lambda, weight=weight, BIC = bic.score, df = df)
  return(out)
}
