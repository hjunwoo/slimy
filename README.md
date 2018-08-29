# slimy
(SampLIng-based Method for causalitY)

## Institute for Health Informatics, University of Minnesota

- slimy is a collection of Markov chain Monte Carlo (MCMC) 
samplers, implemented in R, designed to infer directed acyclic graphs 
representing causal relationships between nodes.
- Both discrete (multinomial) and continuous (Gaussian) data can be treated.
- Default algorithms include Gibbs sampling using parent set-based blocking 
  (Goudie and Mukherjee), which enables efficient sampling of 
  directed acyclic graph (DAG) space for a large number of variables.
