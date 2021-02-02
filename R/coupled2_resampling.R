#' @rdname coupled2_resampling
#' @title (2-way) Coupled multinomial resampling using common uniforms random variables 
#' @description This function performs (2-way) coupled multinomial resampling 
#' on each particle system using common uniform random variables.
#' @param normweights1 normalized weights of first particle system
#' @param normweights2 normalized weights of second particle system
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return a matrix of ancestor indexes of size ndraws x 2.
#' @export
coupled2_resampling <- function(normweights1, normweights2, ndraws, rand){
  ancestors <- matrix(0, nrow = ndraws, ncol = 2)
  ancestors[, 1] <- multinomial_resampling(normweights1, ndraws, rand)
  ancestors[, 2] <- multinomial_resampling(normweights2, ndraws, rand)
  return(ancestors)
  
}

#' @rdname coupled2_maximal_independent_residuals
#' @title (2-way) Maximally coupled multinomial resampling with independent residuals
#' @description This function performs (2-way) maximally coupled multinomial resampling 
#' with independent residuals on each particle system.
#' @param normweights1 normalized weights of first particle system
#' @param normweights2 normalized weights of second particle system
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return a matrix of ancestor indexes of size ndraws x 2.
#' @export
coupled2_maximal_independent_residuals <- function(normweights1, normweights2, ndraws, rand){
  # compute overlap 
  nu <- pmin(normweights1, normweights2)
  alpha <- sum(nu)
  mu <- nu / alpha
  
  # compute residuals
  R1 <- (normweights1 - nu) / (1 - alpha)
  R2 <- (normweights2 - nu) / (1 - alpha)  
  
  # preallocate
  ancestors <- matrix(0, nrow = ndraws, ncol = 2)
  
  # sample from overlap
  uniforms <- runif(ndraws)
  coupled <- (uniforms < alpha)
  ncoupled <- sum(coupled)
  if (ncoupled > 0){
    ancestors[coupled] <- multinomial_resampling(mu, ncoupled, runif(ncoupled))
  }
  
  # sample from residuals independently
  if (ncoupled < ndraws){
    more_uniforms <- runif(ndraws - ncoupled)
    ancestors[!coupled, 1] <- multinomial_resampling(R1, ndraws - ncoupled, more_uniforms)
    
    more_uniforms <- runif(ndraws - ncoupled)
    ancestors[!coupled, 2] <- multinomial_resampling(R2, ndraws - ncoupled, more_uniforms)
  }
  
  return(ancestors)
  
}

#' @rdname coupled2_maximal_coupled_residuals
#' @title (2-way) Maximally coupled multinomial resampling with coupled residuals 
#' using common uniform random variables
#' @description This function performs (2-way) maximally coupled multinomial resampling 
#' with coupled residuals using common uniform random variables on each particle system.
#' @param normweights1 normalized weights of first particle system
#' @param normweights2 normalized weights of second particle system
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return a matrix of ancestor indexes of size ndraws x 2.
#' @export
coupled2_maximal_coupled_residuals <- function(normweights1, normweights2, ndraws, rand){
  # compute overlap 
  nu <- pmin(normweights1, normweights2)
  alpha <- sum(nu)
  mu <- nu / alpha
  
  # compute residuals
  R1 <- (normweights1 - nu) / (1 - alpha)
  R2 <- (normweights2 - nu) / (1 - alpha)  
  
  # preallocate
  ancestors <- matrix(0, nrow = ndraws, ncol = 2)
  
  # sample from overlap
  uniforms <- runif(ndraws)
  coupled <- (uniforms < alpha)
  ncoupled <- sum(coupled)
  if (ncoupled > 0){
    ancestors[coupled] <- multinomial_resampling(mu, ncoupled, runif(ncoupled))
  }
  
  # sample from residuals using common uniform random variables
  more_uniforms <- runif(ndraws - ncoupled)
  if (ncoupled < ndraws){
    ancestors[!coupled] <- coupled2_resampling(R1, R2, ndraws - ncoupled, more_uniforms) 
  }
  
  return(ancestors)
  
}