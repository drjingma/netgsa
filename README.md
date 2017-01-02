# netgsa
Network-based Gene Set Analysis

1. The current version gives an error message if either network has no interactions at all, but this should be allowed in future versions. 

2. We should also include one sample test using NetGSA.

3. For some settings, the variance estimates from approxVarEst could be negative. This affects the subsequent estimation with profile likelihood. We should thus add a validity check after approxVarEst. If negative, then use the previous simple estimates as the initialization for profileVarEst. In addition, approxVarEst may not converge for p > 2000. In such settings, we directly apply profileVarEst. 
