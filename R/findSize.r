#' Determine a sample size
#'
#' Determine a sample size for a given power and hyperparameter configuration
#' 
#' @param power desired power
#' @param a1 alpha, the hyperparameter of the prior distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the prior distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the prior distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the prior distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the prior distribution under the null
#' @param b beta, the hyperparameter of the prior distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param family "binomial" or "poisson", depending on test
#' @param method "ones" or "halves"
#' @seealso \code{\link{samplePower}}
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' 	
#' a1 <- 2; b1 <- 6
#' a2 <- 3; b2 <- 3
#' findSize(power = .8, a1, b1, a2, b2, family = "poisson")
#' 
#' a1 <- 5; b1 <- .2 # E = 1
#' a2 <- 30; b2 <- .167 # E = 5
#' findSize(power = .8, a1, b1, a2, b2, family = "poisson")
#' 
#' a1 <- 4; b1 <- 2
#' a2 <- 6; b2 <- 2
#' findSize(power = .5, a1, b1, a2, b2)
#' findSize(power = .6, a1, b1, a2, b2)
#' findSize(power = .7, a1, b1, a2, b2)
#' findSize(power = .8, a1, b1, a2, b2) # takes a while
#' findSize(power = .8, a1, b1, a2, b2, method = "halves")
#' 
#' findSize(power = .9, a1, b1, a2, b2, method = "halves")
#' 
#' 
#' }
#'
findSize <- function(power, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1, 
  family = c("binomial", "poisson"),
  method = c("ones", "halves")
){
	
  # check arguments
  family <- match.arg(family)  
  method <- match.arg(method)  
  
  # make a convenience function
  last <- function(v) tail(v, 1)
	
  # initialize
  count <- 0
  t <- 2^count
  p <- samplePower(t, a1, b1, a2, b2, a = a, b = b, 
      pi0 = pi0, pi1 = pi1, c = c, family = family)
  if(p > power){
    message(paste0(
      "size = ", 1, 
      ", power = ", round(last(p), 4))
    )
    return(list(size = 1, power = last(p)))
  }
  
  # start routine, march forward  
  while(last(p) < power){
    count <- count + 1    
    t <- c(t, 2^count)
    p <- c(p, 
      samplePower(2^count, a1, b1, a2, b2, a = a, b = b, 
        pi0 = pi0, pi1 = pi1, c = c, family = family
      )
    )
    message(paste0(
      "size = ", 2^count, 
      ", power = ", round(last(p), 4))
    )
  }
  
  
  if(method == "ones"){
    for(k in (2^count-1):(2^(count-1)+1)){
      candidate_p <- samplePower(k, a1, b1, a2, b2, a = a, b = b, 
        pi0 = pi0, pi1 = pi1, c = c, family = family
      )
      message(paste0(
        "size = ", k, 
        ", power = ", round(candidate_p, 4))
      )      
      if(candidate_p < power){
        return(list(size = k+1, power = last(p)))
      } else {
        t <- c(t, k)
        p <- c(p, candidate_p)
      }      
    }
    return(list(size = k+1, power = last(p)))
  }
  
  

  
  if(method == "halves"){
  
    next_t <- 2^(count-1) + 2^(count-2)
    t <- c(t, next_t)
  
    for(k in 2:(count-1)){
      t_p <- samplePower(next_t, a1, b1, a2, b2, a = a, b = b, 
        pi0 = pi0, pi1 = pi1, c = c, family = family
      )
      message(paste0(
        "size = ", next_t, 
        ", power = ", round(t_p, 4))
      )      
      
      if(t_p < power){
        next_t <- next_t + 2^(count-(k+1))
      } else {
        next_t <- next_t - 2^(count-(k+1))
      }    
      
      t <- c(t, next_t)
      p <- c(p, t_p)
      
    }
    
    if(t_p < power){
      #next_t <- next_t+1
      t_p <- samplePower(next_t, a1, b1, a2, b2, a = a, b = b, 
        pi0 = pi0, pi1 = pi1, c = c, family = family
      )
      message(paste0(
        "size = ", next_t, 
        ", power = ", round(t_p, 4))
      )  
      t <- c(t, next_t)
      p <- c(p, t_p)
    }
    
    return(list(size = last(t), power = last(p)))
  }
  
  
  
  
  
}





