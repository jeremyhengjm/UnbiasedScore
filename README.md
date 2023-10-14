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

## Getting started
The two key functions which implements randomized multilevel Monte Carlo and return unbiased estimators of [Rhee and Glynn](https://doi.org/10.1287/opre.2015.1404) are 
`independent_sum` (as considered in the article) and `single_term` (a simpler alternative not considered in the article). These two functions require the following input arguments:
- `model`: a list of objects describing the hidden Markov model obtained by time-discretizing the diffusion process (see below for detailed description);
- `theta`: a vector of model parameters;
- `observations`: a matrix of observations, with rows given by the number of observations, and columns given by the dimension of the observations;
- `nparticles`: an integer specifying number of particles in the conditional particle filter (CPF);
- `resampling_threshold`: a numeric value between zero to one defining the effective sample size threshold below which resampling is triggered at each observation time;
- `coupled2_resampling`: a function defining a coupled resampling scheme for 2 CPFs, setting as `coupled2_maximal_independent_residuals` employs Algorithm 3 of the article;
- `coupled4_resampling`: a function defining a coupled resampling scheme for 4 CPFs, setting as `coupled4_maximalchains_maximallevels_independent_residuals` employs Algorithm 5 of the article;
- `initialization`: a character specifying choice of distribution to initialize chain, taken as `dynamics` for the law of the time-discretized diffusion dynamics 
or `particlefilter` for the law of a trajectory sampled from a particle filter;
- `algorithm`: a character specifying the type of conditional particle filter, taken as `CPF` for the standard CPF, `CASPF` for CPF with ancestor sampling, and `CBSPF` for CPF with backward sampling;
- `k`: an integer specifying a tuning parameter of our unbiased estimator which controls the iteration at which to start averaging, taken as zero by default;
- `m`: an integer specifying a tuning parameter of our unbiased estimator which controls the iteration at which to stop averaging, taken as one by default;
- `level_distribution`: a list containing `mass_function` and `tail_function` that specify the distribution of levels, e.g. by calling the function `compute_level_distribution`; 

The `model` input is given by the output of the function `hmm_ornstein_uhlenbeck` for Example 1, the function `hmm_logistic_diffusion_full` for Example 2, 
and the function `hmm_neuroscience_diffusion` for Example 3. Users who want to apply this package to other problems will have to write a function defining 
the hidden Markov model that returns the following objects:
- `xdimension`: an integer specifying the dimension of the latent diffusion process;
- `ydimension`: an integer specifying the dimension of the observation process;
- `theta_dimension`: an integer specifying the dimension of the parameter space;
- `theta_names`: a vector of characters to enumerate the parameters;
- `theta_positivity`: a vector of logicals to index parameters with positivity constraints;
- `construct_discretization`: a function constructing time-discretization objects for a given discretization level (see the implementation in existing functions);
- `construct_successive_discretization`: a function constructing time-discretization objects for two successive discretization levels (see the implementation in existing functions);
- `sigma`: an object defining the diffusion coefficient of the process; 
- `rinit`: a function to sample from the initial distribution of the diffusion process; 
- `rtransition`: a function to sample from the Markov transition defined by the time-discretized diffusion process;
- `dtransition`: a function to evaluate the the Markov transition defined by the time-discretized diffusion process;
- `dmeasurement`: a function to evaluate the observation/measurement density; 
- `functional`: a function whose conditional expectation one wants to estimate unbiasedly.

To get familarized with the package and its numerical examples, users are strongly encouraged to go through the R notebooks `ornstein_uhlenbeck.Rmd`, `logistic_diffusion.Rmd`, and `neuroscience_diffusion.Rmd` in the folder `notebooks/`. The demo on the Ornstein-Uhlenbeck model should be read first as it serves to guide the reader through various functionalities of the package on a simple example. 



