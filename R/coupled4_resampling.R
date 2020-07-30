#' @rdname coupled4_resampling
#' @title (4-way) Coupled multinomial resampling using common uniforms random variables 
#' @description This function performs (4-way) coupled multinomial resampling 
#' on each particle system using common uniform random variables
#' @param normweights1 normalized weights of first particle system
#' @param normweights2 normalized weights of second particle system
#' @param normweights3 normalized weights of third particle system
#' @param normweights4 normalized weights of fourth particle system
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4
#' @export
coupled4_resampling <- function(normweights1, normweights2, normweights3, normweights4, ndraws, rand){
  ancestors <- matrix(0, nrow = ndraws, ncol = 4)
  ancestors[, 1] <- multinomial_resampling(normweights1, ndraws, rand)
  ancestors[, 2] <- multinomial_resampling(normweights2, ndraws, rand)
  ancestors[, 3] <- multinomial_resampling(normweights3, ndraws, rand)
  ancestors[, 4] <- multinomial_resampling(normweights4, ndraws, rand)
  return(ancestors)
  
}

#' @rdname coupled4_maximal_independent_residuals
#' @title (4-way) Maximally coupled multinomial resampling with independent residuals
#' @description This function performs (4-way) maximally coupled multinomial resampling 
#' with independent residuals on each particle system
#' @param normweights1 normalized weights of first particle system
#' @param normweights2 normalized weights of second particle system
#' @param normweights3 normalized weights of third particle system
#' @param normweights4 normalized weights of fourth particle system
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4s
#' @export
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

#' @rdname coupled4_maximal_coupled_residuals
#' @title (4-way) Maximally coupled multinomial resampling with coupled residuals 
#' using common uniform random variables
#' @description This function performs (4-way) maximally coupled multinomial resampling 
#' with coupled residuals using common uniform random variables on each particle system
#' @param normweights1 normalized weights of first particle system
#' @param normweights2 normalized weights of second particle system
#' @param normweights3 normalized weights of third particle system
#' @param normweights4 normalized weights of fourth particle system
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4
#' @export
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

#' @rdname coupled4_maximallevels_maximalchains_independent_residuals
#' @title (4-way) coupled multinomial resampling by maximally coupling pairs of indexes 
#' across levels and maximally coupling chains on each level with independent residuals 
#' @description This function performs (4-way) coupled multinomial resampling by maximally coupling pairs 
#' of indexes across levels and maximally coupling chains on each level with independent residuals 
#' on each particle system
#' @param normweights_coarse1 normalized weights of first particle system on coarse level
#' @param normweights_coarse2 normalized weights of second particle system on coarse level
#' @param normweights_fine1 normalized weights of first particle system on fine level
#' @param normweights_fine2 normalized weights of second particle system on fine level
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4
#' @export
coupled4_maximallevels_maximalchains_independent_residuals <- function(normweights_coarse1, normweights_coarse2, 
                                                                       normweights_fine1, normweights_fine2, 
                                                                       ndraws, rand){
  
  output <- maximal_maximal_multinomial_resampling(normweights_coarse1, normweights_coarse2, 
                                                   normweights_fine1, normweights_fine2, 
                                                   ndraws, residualtype = 1)
  
  ancestors <- matrix(nrow = ndraws, ncol = 4)
  ancestors[, 1] <- output[, 1]
  ancestors[, 2] <- output[, 2]
  ancestors[, 3] <- output[, 3]
  ancestors[, 4] <- output[, 4]
  
  return(ancestors)
  
}

#' @rdname coupled4_maximallevels_maximalchains_coupled_residuals
#' @title (4-way) coupled multinomial resampling by maximally coupling pairs of indexes 
#' across levels and maximally coupling chains on each level with 
#' coupled residuals using common uniform random variables 
#' @description This function performs (4-way) coupled multinomial resampling by maximally coupling pairs 
#' of indexes across levels and maximally coupling chains on each level with with coupled residuals 
#' using common uniform random variables on each particle system 
#' @param normweights_coarse1 normalized weights of first particle system on coarse level
#' @param normweights_coarse2 normalized weights of second particle system on coarse level
#' @param normweights_fine1 normalized weights of first particle system on fine level
#' @param normweights_fine2 normalized weights of second particle system on fine level
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4
#' @export
coupled4_maximallevels_maximalchains_coupled_residuals <- function(normweights_coarse1, normweights_coarse2, 
                                                                   normweights_fine1, normweights_fine2, 
                                                                   ndraws, rand){
  
  output <- maximal_maximal_multinomial_resampling(normweights_coarse1, normweights_coarse2, 
                                                   normweights_fine1, normweights_fine2, 
                                                   ndraws, residualtype = 2)
  
  ancestors <- matrix(nrow = ndraws, ncol = 4)
  ancestors[, 1] <- output[, 1]
  ancestors[, 2] <- output[, 2]
  ancestors[, 3] <- output[, 3]
  ancestors[, 4] <- output[, 4]
  
  return(ancestors)
  
}

#' @rdname coupled4_maximalchains_maximallevels_independent_residuals
#' @title (4-way) coupled multinomial resampling by maximally coupling pairs of chains on each level 
#' and maximally coupling across levels with independent residuals 
#' @description This function performs (4-way) coupled multinomial resampling by maximally coupling pairs of chains  
#' and maximally coupling across levels with independent residuals on each particle system 
#' @param normweights_coarse1 normalized weights of first particle system on coarse level
#' @param normweights_coarse2 normalized weights of second particle system on coarse level
#' @param normweights_fine1 normalized weights of first particle system on fine level
#' @param normweights_fine2 normalized weights of second particle system on fine level
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4
#' @export
coupled4_maximalchains_maximallevels_independent_residuals <- function(normweights_coarse1, normweights_coarse2, 
                                                                       normweights_fine1, normweights_fine2, 
                                                                       ndraws, rand){
  
  ancestors <- matrix(nrow = ndraws, ncol = 4)
  
  if (all(normweights_coarse1 == normweights_coarse2) & !all(normweights_fine1 == normweights_fine2)) {
    # chains on coarse level have met and chains on fine have not met
    # sample maximal coupling between levels for first pair
    output <- coupled2_maximal_independent_residuals(normweights_coarse1, normweights_fine1, ndraws, rand)
    ancestors[, 1] <- output[, 1]
    ancestors[, 3] <- output[, 2]
    
    # assign the same ancestor index on the coarse level
    ancestors[, 2] <- ancestors[, 1]
    
    # rejection sampling to sample from maximal coupling for the second pair
    ancestors[, 4] <- maximal_rejection_sampling(ancestors[, 2], normweights_coarse2, normweights_fine2)
    
  } else if (!all(normweights_coarse1 == normweights_coarse2) & all(normweights_fine1 == normweights_fine2)) {
    # chains on coarse level have not met and chains on fine level have met
    # sample maximal coupling between levels for first pair
    output <- coupled2_maximal_independent_residuals(normweights_coarse1, normweights_fine1, ndraws, rand)
    ancestors[, 1] <- output[, 1]
    ancestors[, 3] <- output[, 2]
    
    # assign the same ancestor index on the fine level
    ancestors[, 4] <- ancestors[, 3]
    
    # rejection sampling to sample from maximal coupling for the second pair
    ancestors[, 2] <- maximal_rejection_sampling(ancestors[, 4], normweights_fine2, normweights_coarse2)
    
  } else {
    # either no meeting has occurred or chains on both levels have met (faithfulness holds in this case)
    output <- maximal_maximal_multinomial_resampling(normweights_coarse1, normweights_fine1, 
                                                     normweights_coarse2, normweights_fine2, 
                                                     ndraws, residualtype = 1)
    ancestors[, 1] <- output[, 1]
    ancestors[, 2] <- output[, 3]
    ancestors[, 3] <- output[, 2]
    ancestors[, 4] <- output[, 4]
  }
  
  return(ancestors)
}

#' @rdname coupled4_maximalchains_maximallevels_coupled_residuals
#' @title (4-way) coupled multinomial resampling by maximally coupling pairs of chains on each level 
#' and maximally coupling across levels with coupled residuals using common uniform random variables 
#' @description This function performs (4-way) coupled multinomial resampling by maximally coupling pairs of chains  
#' and maximally coupling across levels with coupled residuals using common uniform random variables on each particle system
#' @param normweights_coarse1 normalized weights of first particle system on coarse level
#' @param normweights_coarse2 normalized weights of second particle system on coarse level
#' @param normweights_fine1 normalized weights of first particle system on fine level
#' @param normweights_fine2 normalized weights of second particle system on fine level
#' @param ndraws number of samples from multinomial distribution
#' @param rand common uniform random variables used for the coupled resampling
#' @return A matrix of ancestor indexes of size ndraws x 4
#' @export
coupled4_maximalchains_maximallevels_coupled_residuals <- function(normweights_coarse1, normweights_coarse2, 
                                                                   normweights_fine1, normweights_fine2, 
                                                                   ndraws, rand){
  
  ancestors <- matrix(nrow = ndraws, ncol = 4)
  
  if (all(normweights_coarse1 == normweights_coarse2) & !all(normweights_fine1 == normweights_fine2)) {
    # chains on coarse level have met and chains on fine have not met
    # sample maximal coupling between levels for first pair
    output <- coupled2_maximal_coupled_residuals(normweights_coarse1, normweights_fine1, ndraws, rand)
    ancestors[, 1] <- output[, 1]
    ancestors[, 3] <- output[, 2]
    
    # assign the same ancestor index on the coarse level
    ancestors[, 2] <- ancestors[, 1]
    
    # rejection sampling to sample from maximal coupling for the second pair
    ancestors[, 4] <- maximal_rejection_sampling(ancestors[, 2], normweights_coarse2, normweights_fine2)
    
  } else if (!all(normweights_coarse1 == normweights_coarse2) & all(normweights_fine1 == normweights_fine2)) {
    # chains on coarse level have not met and chains on fine level have met
    # sample maximal coupling between levels for first pair
    output <- coupled2_maximal_coupled_residuals(normweights_coarse1, normweights_fine1, ndraws, rand)
    ancestors[, 1] <- output[, 1]
    ancestors[, 3] <- output[, 2]
    
    # assign the same ancestor index on the fine level
    ancestors[, 4] <- ancestors[, 3]
    
    # rejection sampling to sample from maximal coupling for the second pair
    ancestors[, 2] <- maximal_rejection_sampling(ancestors[, 4], normweights_fine2, normweights_coarse2)
    
  } else {
    # either no meeting has occurred or chains on both levels have met (faithfulness holds in this case)
    output <- maximal_maximal_multinomial_resampling(normweights_coarse1, normweights_fine1, 
                                                     normweights_coarse2, normweights_fine2, 
                                                     ndraws, residualtype = 2)
    ancestors[, 1] <- output[, 1]
    ancestors[, 2] <- output[, 3]
    ancestors[, 3] <- output[, 2]
    ancestors[, 4] <- output[, 4]
  }
  
  return(ancestors)
}