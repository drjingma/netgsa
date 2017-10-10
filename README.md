# netgsa
Network-based Gene Set Analysis

1. The current version gives an error message if either network has no interactions at all, but the case where only one network is empty should be allowed in future versions. (If one network is empty, it reduces to a fixed effects model. There's no need to apply this network-based method. Report an error.)

2. We should also include one sample test using NetGSA.

3. For some settings, the variance estimates from approxVarEst could be negative. This affects the subsequent estimation with profile likelihood. We should thus add a validity check after approxVarEst. If negative, then use the previous simple estimates as the initialization for profileVarEst. In addition, approxVarEst may not converge for p > 2000. In such settings, we directly apply profileVarEst. Ali also mentioned potential identifiability issue with NetGSA. (Hopefully the new proposal provides a better initial estimate for use in NetGSA.)

4. When do we get an empty estimated network based on a pre-specified graph and data? What guidance to provide in such a setting? (Return some warning message if this happens. This can't be used in NetGSA.)

5. New error message when running permuted NetGSA: 
    
    Error in bic.netEst(X = t(current_data[[k]]), one = oneMat, lambda = lambda_vec[i],  :
      task 1 failed - "'to' must be of length 1"
    Calls: %dopar% -> 
    
    If we look further at the details, the error happens because:
    
    Error in seq.default(1, dim(zero.pos)[1]) : 'to' must be of length 1
 
    This happens when the input zero matrix is empty. I've made changes to the offline file.
6. That the input matrices for NetGSA and network estimation are of different dimensions, one is p by n, whereas the other is n by p. 

7. Need to make it clear that NetGSA is designed for inference with large sample sizes. If the sample size is too small, return a warning message. 
    
