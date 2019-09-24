adj2inf <-
function(AA, isDAG = FALSE, eta=0.1){
    p = dim(AA)[1]
    Ip = diag( rep(1,p) )
    
    if(isDAG){
      
      DD = solveCpp(Ip - AA)
      return(DD)
      
    }else{
      tmp = eigen(AA, symmetric=FALSE, only.values=TRUE)$values
      
      if( max(abs(tmp)) < 1-eta ){
        DD = solveCpp(Ip - AA)
      }else{
        # This is consistent with the ref #2
        tmp = apply( abs(AA), 1, sum )
        tmp = matrix(tmp+eta, ncol=p, nrow=p, byrow=FALSE) 
        tmp = AA / tmp
        DD = solveCpp(Ip - tmp)
      }  
      
      return(DD)
    }
  }
