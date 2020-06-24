newton <-
function(x0, lb, ub, f, g, h, alpha = 0.25, beta = 0.5, max.iter = 100, tol = 1e-2){
    count = 0  
    x = x0
    
    repeat {
      count = count + 1
      
      ## newton's step
      delta = - g(x)/h(x)
      
      ## line search to pick the stepsize 
      size = 1
      lhs = f(x + size * delta)
      rhs = f(x) + alpha * size * g(x) * delta
      while ( (x + size * delta <=0) || lhs > rhs || (log(lhs) > log(rhs)) ) {
        size = beta * size
        #Trying to minimize evaluation of f() and g():
        if ( !(x + size * delta <=0) ){
          lhs = f(x + size * delta)
          rhs = f(x) + alpha * size * g(x) * delta
        }
      }
      ## Update
      x.new = x + size * delta
      
      if ( count >= max.iter || abs(x-x.new)< tol || x.new > ub || x.new < lb ){
        
        if(count == max.iter) warning("Maximum number of iterations reached!")
        break
      }
      x = x.new
      
    }  
    
    return(list(solution = x, iter = count, stepSize = size))
  }
