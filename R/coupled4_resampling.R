#'@rdname coupled4_resampling
#'@title (4-way) Coupled multinomial resampling using common uniforms random variables 
#'@description This function performs (4-way) coupled multinomial resampling 
#'on each particle system using common uniform random variables
#'@param normweights1 normalized weights of first particle system
#'@param normweights2 normalized weights of second particle system
#'@param normweights3 normalized weights of third particle system
#'@param normweights4 normalized weights of fourth particle system
#'@param ndraws number of samples from multinomial distribution
#'@param rand common uniform random variables used for the coupled resampling
#'@return A matrix of ancestor indexes of size ndraws x 4
#'@export
coupled2_resampling <- function(normweights1, normweights2, normweights3, normweights4, ndraws, rand){
  ancestors <- matrix(0, nrow = ndraws, ncol = 4)
  ancestors[, 1] <- multinomial_resampling(nornormweights1, ndraws, rand)
  ancestors[, 2] <- multinomial_resampling(nornormweights2, ndraws, rand)
  ancestors[, 3] <- multinomial_resampling(nornormweights3, ndraws, rand)
  ancestors[, 4] <- multinomial_resampling(nornormweights4, ndraws, rand)
  return(ancestors)
  
}

#'@rdname coupled4_maximal_independent_residuals
#'@title (4-way) Maximally coupled multinomial resampling with independent residuals
#'@description This function performs (4-way) maximally coupled multinomial resampling 
#'with independent residuals on each particle system
#'@param normweights1 normalized weights of first particle system
#'@param normweights2 normalized weights of second particle system
#'@param normweights3 normalized weights of third particle system
#'@param normweights4 normalized weights of fourth particle system
#'@param ndraws number of samples from multinomial distribution
#'@param rand common uniform random variables used for the coupled resampling
#'@return A matrix of ancestor indexes of size ndraws x 4s
#'@export
coupled4_maximal_independent_residuals <- function(normweights1, normweights2, normweights3, normweights4, ndraws, rand){
  # compute overlap 
  nu <- pmin(normweights1, normweights2, normweights3, normweights4)
  alpha <- sum(nu)
  mu <- nu / alpha
  
  # compute residuals
  R1 <- (normweights1 - nu) / (1 - alpha)
  R2 <- (normweights2 - nu) / (1 - alpha)  
  R3 <- (normweights3 - nu) / (1 - alpha)  
  R4 <- (normweights4 - nu) / (1 - alpha)  
  
  # preallocate
  ancestors <- matrix(0, nrow = ndraws, ncol = 4)
  
  # sample from overlap
  uniforms <- runif(ndraws)
  coupled <- (uniforms < alpha)
  ncoupled <- sum(coupled)
  if (ncoupled > 0){
    ancestors[coupled] <- multinomial_resampling(mu, ncoupled, runif(ncoupled))
  }
  
  # sample from residuals independently
  if (ncoupled < ndraws){
    
    if (all(normweights1 == normweights2)){
      # check if pairs are equal to ensure faithfulness
      more_uniforms <- runif(ndraws - ncoupled)
      ancestors[!coupled, 1] <- multinomial_resampling(R1, ndraws - ncoupled, more_uniforms)
      ancestors[!coupled, 2] <- ancestors[!coupled, 1]
      
    } else {
      # otherwise sample pairs independently
      more_uniforms <- runif(ndraws - ncoupled)
      ancestors[!coupled, 1] <- multinomial_resampling(R1, ndraws - ncoupled, more_uniforms)
      
      more_uniforms <- runif(ndraws - ncoupled)
      ancestors[!coupled, 2] <- multinomial_resampling(R2, ndraws - ncoupled, more_uniforms)
    }
    
    if (all(normweights3 == normweights4)){
      # check if pairs are equal to ensure faithfulness
      more_uniforms <- runif(ndraws - ncoupled)
      ancestors[!coupled, 3] <- multinomial_resampling(R3, ndraws - ncoupled, more_uniforms)
      ancestors[!coupled, 4] <- ancestors[!coupled, 3]
      
    } else {
      # otherwise sample pairs independently
      more_uniforms <- runif(ndraws - ncoupled)
      ancestors[!coupled, 3] <- multinomial_resampling(R3, ndraws - ncoupled, more_uniforms)
      
      more_uniforms <- runif(ndraws - ncoupled)
      ancestors[!coupled, 4] <- multinomial_resampling(R4, ndraws - ncoupled, more_uniforms)
    }
  }
  
  return(ancestors)
  
}

#'@rdname coupled4_maximal_coupled_residuals
#'@title (4-way) Maximally coupled multinomial resampling with coupled residuals 
#'using common uniform random variables
#'@description This function performs (4-way) maximally coupled multinomial resampling 
#'with coupled residuals using common uniform random variables on each particle system
#'@param normweights1 normalized weights of first particle system
#'@param normweights2 normalized weights of second particle system
#'@param normweights3 normalized weights of third particle system
#'@param normweights4 normalized weights of fourth particle system
#'@param ndraws number of samples from multinomial distribution
#'@param rand common uniform random variables used for the coupled resampling
#'@return A matrix of ancestor indexes of size ndraws x 4
#'@export
coupled4_maximal_coupled_residuals <- function(normweights1, normweights2, normweights3, normweights4, ndraws, rand){
  # compute overlap 
  nu <- pmin(normweights1, normweights2, normweights3, normweights4)
  alpha <- sum(nu)
  mu <- nu / alpha
  
  # compute residuals
  R1 <- (normweights1 - nu) / (1 - alpha)
  R2 <- (normweights2 - nu) / (1 - alpha)  
  R3 <- (normweights3 - nu) / (1 - alpha)  
  R4 <- (normweights4 - nu) / (1 - alpha)  
  
  # preallocate
  ancestors <- matrix(0, nrow = ndraws, ncol = 4)
  
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
    ancestors[!coupled] <- coupled4_resampling(R1, R2, R3, R4, ndraws - ncoupled, more_uniforms) 
  }
  
  return(ancestors)
  
}