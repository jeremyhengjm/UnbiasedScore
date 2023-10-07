# UnbiasedScore
This R package provides unbiased estimators for a class of partially observed diffusions. It facilitates parameter inference using stochastic gradient algorithms 
by considering unbiased estimation of the score function, i.e. the gradient of the log-likelihood function for these models. The implementation is based on the arXiv preprint:

[On Unbiased Score Estimation for Partially Observed Diffusions](https://arxiv.org/abs/2105.04912) by Jeremy Heng, Jeremie Houssineau, and Ajay Jasra. 

## Numerical examples
The methodology is illustrated on three examples:
1. Ornstein--Uhlenbeck process;
2. Logistic diffusion model for population dynamics of red kangaroos;
3. Neural network model for grid cells in the medial entorhinal cortex.

The scripts to reproduce the results in the article are in the `inst` folder, with subfolders corresponding to each example. Run scripts with filenames starting with "simulate" to generate results that are saved in the `inst/example/results` folder. Run scripts with filenames starting with "plot" to create plots of the results. 

Example 1 is considered as it is a simple model where it is possible to compute the score function exactly using Kalman smoothing. 
Example 2 concerns a dataset from [Caughley et al. (1987)](https://doi.org/10.2307/4946) that is stored in `inst/logistic_diffusion/kangaroo.RData`.
Example 3 is an application on the data collected by [Hafting et al. (2008)](https://doi.org/10.1038/nature06957) and is stored in `inst/neuroscience_diffusion/gridcells.RData`.

## R Notebook tutorials


