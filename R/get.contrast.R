get.contrast <-
function(InfMat, b){
  ##InfMat: a list of ncond information matrices
  ##b: the indicator for one pathway (a row vector of 0 and 1's)
  ncond = length(InfMat)
  p = nrow(InfMat[[1]])
  LC = matrix(0, nrow=ncond - 1, ncol = ncond * p)
  for (j in 1:(ncond-1)){
    L1 = (b %*% InfMat[[j]]) * b
    L2 = (b %*% InfMat[[j+1]]) * b
    tmp = c(-L1, L2)
    LC[j, ((j-1)*p + 1): ((j+1)*p)] = tmp 
  }
  return(LC)
}
