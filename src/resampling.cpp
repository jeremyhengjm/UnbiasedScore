#include <RcppEigen.h>

using namespace Rcpp;
using namespace std;  

// [[Rcpp::export]]
IntegerVector multinomial_resampling(const NumericVector & weights, 
                                     int ndraws, 
                                     const NumericVector & rand){
  IntegerVector ancestors(ndraws);
  NumericVector cumsumw = cumsum(weights);
  int i;
  for (int n = 0; n < ndraws; n++){
    i = 0;
    while (cumsumw(i) < rand(n)){
      i = i + 1;
    }
    ancestors(n) = i + 1;
  }
  return ancestors;
}

// [[Rcpp::export]]
IntegerMatrix maximal_maximal_multinomial_resampling(const NumericVector & weights1,
                                                     const NumericVector & weights2,
                                                     const NumericVector & weights3,
                                                     const NumericVector & weights4,
                                                     int ndraws,
                                                     int residualtype){
  
  IntegerMatrix ancestors(ndraws, 4); // output
  int nparticles = weights1.size();
  NumericVector overlap_first(nparticles);
  NumericVector overlap_second(nparticles);
  double size_first = 0;
  double size_second = 0;
  double minvalue;
  int index1, index2, index3, index4;
  double uniform, loguniform;
  IntegerVector output_resampling;
  NumericVector random(1);
  double pmf_first, pmf_second;
  double cumsum_upper1, cumsum_upper2, cumsum_lower1, cumsum_lower2;
  
  NumericVector prob_overlap_first(nparticles);
  NumericVector prob_overlap_second(nparticles);
  NumericVector residual1(nparticles);
  NumericVector residual2(nparticles);
  NumericVector residual3(nparticles);
  NumericVector residual4(nparticles);
  
  double logacceptprob, logrejectprob;
  bool reject;
  // double nattempts;
  
  // Compute overlap and its size
  for (int n = 0; n < nparticles; n++){
    // First pair
    minvalue = min(weights1(n), weights2(n));
    overlap_first(n) = minvalue;
    size_first += minvalue;
    
    // Second pair
    minvalue = min(weights3(n), weights4(n));
    overlap_second(n) = minvalue;
    size_second += minvalue;
  }
  
  // Normalized overlap and marginals of residual
  
  for (int n = 0; n < nparticles; n++){
    prob_overlap_first(n) = overlap_first(n) / size_first;
    prob_overlap_second(n) = overlap_second(n) / size_second;
    
    residual1(n) = (weights1(n) - overlap_first(n)) / (1 - size_first);
    residual2(n) = (weights2(n) - overlap_first(n)) / (1 - size_first);
    
    residual3(n) = (weights3(n) - overlap_second(n)) / (1 - size_second);
    residual4(n) = (weights4(n) - overlap_second(n)) / (1 - size_second);
  }
  
  // Compute cumulative sum of all four residuals needed for common uniform coupling
  NumericVector cumsum_residual1 = cumsum(residual1);
  NumericVector cumsum_residual2 = cumsum(residual2);
  NumericVector cumsum_residual3 = cumsum(residual3);
  NumericVector cumsum_residual4 = cumsum(residual4);
  
  for (int n = 0; n < ndraws; n++){
    // Sample from first maximal coupling
    
    uniform = R::runif(0,1); // this returns a double
    if (uniform < size_first){
      random = runif(1); // this returns a vector of size one
      output_resampling = multinomial_resampling(prob_overlap_first, 1, random);
      index1 = output_resampling(0);
      index2 = output_resampling(0);
      
    } else {
      // Sample first and second index from residuals
      if (residualtype == 1){
        // Sample first and second index independently
        random = runif(1);
        output_resampling = multinomial_resampling(residual1, 1, random);
        index1 = output_resampling(0);
        
        random = runif(1);
        output_resampling = multinomial_resampling(residual2, 1, random);
        index2 = output_resampling(0);
        
      }
      
      if (residualtype == 2){
        // Sample first and second index using common uniform random variable
        random = runif(1);
        
        output_resampling = multinomial_resampling(residual1, 1, random);
        index1 = output_resampling(0);
        
        output_resampling = multinomial_resampling(residual2, 1, random);
        index2 = output_resampling(0);
      }
    }
    
    // Evaluate maximal couplings of both pairs at sampled index1 and index2 
    pmf_first = 0;
    pmf_second = 0;
    // Overlaps
    if (index1 == index2){
      pmf_first += overlap_first(index1 - 1); // note that Cpp index begins from 0 
      pmf_second += overlap_second(index1 - 1);
    }
    
    // Residuals
    if (residualtype == 1){
      // Independent coupling of residuals
      pmf_first += (1.0 - size_first) * residual1(index1 - 1) * residual2(index2 - 1);
      pmf_second += (1.0 - size_second) * residual3(index1 - 1) * residual4(index2 - 1);
    }
    
    if (residualtype == 2){
      // Common uniform coupling of residuals
      
      // First pair 
      cumsum_upper1 = cumsum_residual1(index1 - 1);
      cumsum_upper2 = cumsum_residual2(index2 - 1);
      
      if (index1 == 1){
        cumsum_lower1 = 0;
      } else {
        cumsum_lower1 = cumsum_residual1(index1 - 2);
      }
      
      if (index2 == 1){
        cumsum_lower2 = 0;
      } else {
        cumsum_lower2 = cumsum_residual2(index2 - 2);
      }
      
      pmf_first += (1.0 - size_first) * ( min(cumsum_upper1, cumsum_upper2) - max(cumsum_lower1, cumsum_lower2) );
      
      // Second pair 
      cumsum_upper1 = cumsum_residual3(index1 - 1);
      cumsum_upper2 = cumsum_residual4(index2 - 1);
      
      if (index1 == 1){
        cumsum_lower1 = 0;
      } else {
        cumsum_lower1 = cumsum_residual3(index1 - 2);
      }
      
      if (index2 == 1){
        cumsum_lower2 = 0;
      } else {
        cumsum_lower2 = cumsum_residual4(index2 - 2);
      }
      
      pmf_second += (1.0 - size_second) * ( min(cumsum_upper1, cumsum_upper2) - max(cumsum_lower1, cumsum_lower2) );
    }
    
    // Attempt to sample from overlapping pairs
    loguniform = log(R::runif(0,1)); // this returns a double
    logacceptprob = log(pmf_second) - log(pmf_first);
    if (loguniform < logacceptprob){
      index3 = index1;
      index4 = index2;
    } else {
      // Rejection sampler for second pair
      reject = true;
      // nattempts = 0;
      while (reject){
        // nattempts += 1;
        // Rcpp::Rcout << "Number of attempts:" << std::endl << nattempts << std::endl;
        // Sample proposal from second maximal coupling
        uniform = R::runif(0,1); // this returns a double
        if (uniform < size_second){
          // Sample third and fourth index from overlap
          random = runif(1); // this returns a vector of size one
          output_resampling = multinomial_resampling(prob_overlap_second, 1, random);
          index3 = output_resampling(0);
          index4 = output_resampling(0);
          
        } else {
          // Sample third and fourth index from residuals
          if (residualtype == 1){
            // Sample third and fourth index independently
            random = runif(1);
            output_resampling = multinomial_resampling(residual3, 1, random);
            index3 = output_resampling(0);
            
            random = runif(1);
            output_resampling = multinomial_resampling(residual4, 1, random);
            index4 = output_resampling(0);
            
          }
          
          if (residualtype == 2){
            // Sample third and fourth index using common uniform random variable
            random = runif(1);
            
            output_resampling = multinomial_resampling(residual3, 1, random);
            index3 = output_resampling(0);
            
            output_resampling = multinomial_resampling(residual4, 1, random);
            index4 = output_resampling(0);
          }
        }
        
        // Evaluate maximal couplings of both pairs at sampled index3 and index4
        pmf_first = 0;
        pmf_second = 0;
        // Overlaps
        if (index3 == index4){
          pmf_first += overlap_first(index3 - 1); // note that Cpp index begins from 0 
          pmf_second += overlap_second(index3 - 1);
        }
        
        // Residuals
        if (residualtype == 1){
          // Independent coupling of residuals
          pmf_first += (1.0 - size_first) * residual1(index3 - 1) * residual2(index4 - 1);
          pmf_second += (1.0 - size_second) * residual3(index3 - 1) * residual4(index4 - 1);
        }
        
        if (residualtype == 2){
          // Common uniform coupling of residuals
          
          // First pair 
          cumsum_upper1 = cumsum_residual1(index3 - 1);
          cumsum_upper2 = cumsum_residual2(index4 - 1);
          
          if (index3 == 1){
            cumsum_lower1 = 0;
          } else {
            cumsum_lower1 = cumsum_residual1(index3 - 2);
          }
          
          if (index4 == 1){
            cumsum_lower2 = 0;
          } else {
            cumsum_lower2 = cumsum_residual2(index4 - 2);
          }
          
          pmf_first += (1.0 - size_first) * ( min(cumsum_upper1, cumsum_upper2) - max(cumsum_lower1, cumsum_lower2) );
          
          // Second pair 
          cumsum_upper1 = cumsum_residual3(index3 - 1);
          cumsum_upper2 = cumsum_residual4(index4 - 1);
          
          if (index3 == 1){
            cumsum_lower1 = 0;
          } else {
            cumsum_lower1 = cumsum_residual3(index3 - 2);
          }
          
          if (index4 == 1){
            cumsum_lower2 = 0;
          } else {
            cumsum_lower2 = cumsum_residual4(index4 - 2);
          }
          
          pmf_second += (1.0 - size_second) * ( min(cumsum_upper1, cumsum_upper2) - max(cumsum_lower1, cumsum_lower2) );
        }
        
        logrejectprob = log(pmf_first) - log(pmf_second);
        
        // accept or reject
        loguniform = log(R::runif(0,1)); // this returns a double
        if (loguniform < logrejectprob){
          reject = true;
        } else {
          reject = false;
        }
      }
    }
    ancestors(n, 0) = index1; 
    ancestors(n, 1) = index2; 
    ancestors(n, 2) = index3; 
    ancestors(n, 3) = index4; 
  }
  return ancestors;
}