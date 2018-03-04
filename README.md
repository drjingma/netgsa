# netgsa
Network-based Gene Set Analysis

1. There is some potential identifiability issue with NetGSA. For example, the estimates of s2g and s2e can flip.

2. When do we get an empty estimated network based on a pre-specified graph and data? What guidance to provide in such a setting? (Return some warning message if this happens. This can't be used in NetGSA.)

3. We have improved the interface of NetGSA for biomedical scientists. In particular, we have connected netgsa with the bioconductor package KEGGgraph, which allows the input of network information in the form of txt files. In the future we may make NetGSA a built-in option in http://www.metaboanalyst.ca/. The latter requires a bit more work in finding the right KEGG metabolic network, because individual metabolomic studies may not detect all metabolites that are connected in KEGG. Ali has a student working on this. 
    
