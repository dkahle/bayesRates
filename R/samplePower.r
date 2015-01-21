#' Determine the power of a test with given sample size
#'
#' Determine the power of a test with given sample size
#' 
#' @param n the sample size (this is the t parameter in the poisson case)
#' @param a1 alpha, the hyperparameter of the prior distribution
#' @param b1 beta, the hyperparameter of the prior distribution
#' @param a2 alpha, the hyperparameter of the prior distribution
#' @param b2 beta, the hyperparameter of the prior distribution
#' @param a alpha, the hyperparameter of the prior distribution under the null
#' @param b beta, the hyperparameter of the prior distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param family "binomial" or "poisson", depending on test
#' @seealso \code{\link{samplePowerBinomial}}, \code{\link{samplePowerPoisson}}, \code{\link{samplePowerEst}}, \code{\link{find_size}}
#' @export samplePower
#' @examples
#' \dontrun{
#' 
#' a1 <- 3; b1 <- 7
#' a2 <- 7; b2 <- 3
#' plot_beta(c(a1, a2), c(b1, b2))
#' samplePower(n = 10, a1, b1, a2, b2)
#' samplePower(n = 11, a1, b1, a2, b2)
#' samplePower(n = 1, a1, b1, a2, b2)
#' 
#' samplePower(100, a1, b1, a2, b2)
#' samplePowerBinomial(100, a1, b1, a2, b2)
#'
#'
#'
#'
#'
#' # monte carlo check of samplePower
#' N <- 10000
#' n <- 100
#' th1 <- rbeta(N, a1, b1)
#' y1s <- rbinom(N, n, th1)
#' th2 <- rbeta(N, a2, b2)
#' y2s <- rbinom(N, n, th2)
#' a <- a1; b <- b1 # belief under null
#' c <- 1 # relative penalty of type II to type I error
#' pi0 <- .5; pi1 <- .5 # prior likelihoods of hypotheses
#' 
#' samplePower(n = 100, a1, b1, a2, b2)
#' # exact power should be 91.95%
#'
#' test <- function(x) bayes_binom_test(x, n, a1, b1, a2, b2)$reject
#' mean(apply(cbind(y1s, y2s), 1, test))
#' # simulated power checks out at ~92%
#'   
#' 
#' 
#' # check the probability of type I error of the test
#' y2s <- rbinom(N, n, th1) # note : not resampled betas
#' mean(apply(cbind(y1s, y2s), 1, test))
#' # prob(type I error) ~= 2.5% 
#' 
#' # compare the power of the bayesian-derived test to 
#' # that of the standard two-tailed two-sample binomial test for the 
#' # same probability of type I error
#' classic_test <- function(x) prop.test(x, c(n, n))$p.value < .025
#' 
#' # confirm it's probability of type I error
#' mean(apply(cbind(y1s, y2s), 1, classic_test)) # it's low...
#' 
#' # check it's power against the same alternative
#' th2 <- rbeta(N, a2, b2)
#' y2s <- rbinom(N, n, th2)
#' mean(apply(cbind(y1s, y2s), 1, classic_test))
#' # power ~= 87.5%, lower than the bayesian test of 91.52%
#' 
#' 
#' 
#' 
#' 
#' 
#'
#'
#' a1 <- 1; b1 <- 1 # E = 1
#' a2 <- 3; b2 <- 1 # E = 3
#' plot_gamma(c(a1, a2), c(b1, b2))
#' samplePower(n = 2, a1, b1, a2, b2, family = "poisson")
#' samplePower(n = 10, a1, b1, a2, b2, family = "poisson")
#' samplePower(n = 11, a1, b1, a2, b2, family = "poisson")
#' samplePower(n = 10:11, a1, b1, a2, b2, family = "poisson")
#' 
#' samplePower(n = 250, a1, b1, a2, b2, family = "poisson") 
#' samplePowerEst(n = c(250, 500, 1000), a1, b1, a2, b2) # 5s, off by .001
#' 
#' 
#' 
#' a1 <- 4; b1 <- 2
#' a2 <- 6; b2 <- 2
#' plot_gamma(c(a1, a2), c(b1, b2))
#' samplePower(n = 20, a1, b1, a2, b2, family = "poisson")
#' 
#' 
#' }
samplePower <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1,
  family = c("binomial", "poisson")
){

  # check arguments
  stop2 <- function(x) stop(x, call. = FALSE)
  if(missing(n)) stop2("the sample size n must be specified.  see ?samplePower.")
  if(missing(a1)) stop2("the hyperparameter a1 must be specified.  see ?samplePower.")  
  if(missing(b1)) stop2("the hyperparameter b1 must be specified.  see ?samplePower.")  
  if(missing(a2)) stop2("the hyperparameter a2 must be specified.  see ?samplePower.")  
  if(missing(b2)) stop2("the hyperparameter b2 must be specified.  see ?samplePower.")        
  family <- match.arg(family)
  stopifnot( pi0 >= 0 )
  stopifnot( pi1 >= 0 )  
  stopifnot( pi0 + pi1 == 1 )

  # dispatch correct function
  if(family == "binomial"){
    out <- samplePowerBinomial(n, a1, b1, a2, b2, a, b, pi0, pi1, c)
  }

  if(family == "poisson"){
    out <- samplePowerPoisson(n, a1, b1, a2, b2, a, b, pi0, pi1, c)    
  }

  # return
  out
}









































#' Compute the power of the binomial proportions test for a given sample size (n) with specified beta priors
#'
#' Compute the power of the binomial proportions test for a given sample size (n) with specified beta priors
#' 
#' @param n the sample size 
#' @param a1 alpha, the hyperparameter of the beta distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the beta distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the beta distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the beta distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the beta distribution under the null
#' @param b beta, the hyperparameter of the beta distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @seealso \code{\link{samplePower}}, \code{\link{samplePowerEst}}, \code{\link{find_size}}
#' @export samplePowerBinomial
#' @examples
#' \dontrun{
#' 
#' # see ?samplePower
#' 
#' 
#' }
samplePowerBinomial <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1)
{

  # vectorize for multiple n
  if(is.numeric(n) && length(n) > 1){
    return(sapply(as.list(n), function(x){
      message(".", appendLF = FALSE)
      samplePowerBinomial(x, a1, b1, a2, b2, a, b, pi0, pi1, c)
    }))
  }
  
  # Rcpp
  samplePowerBinomialCpp(n, a1, b1, a2, b2, a, b, pi0, pi1, c)
}





















































#' Compute the power of the Poisson rate test for a given sample size (t) with specified gamma priors
#'
#' Compute the power of the Poisson rate test for a given sample size (t) with specified gamma priors
#' 
#' @param t the sample size 
#' @param a1 alpha, the hyperparameter of the gamma distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the gamma distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the gamma distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the gamma distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the gamma distribution under the null
#' @param b beta, the hyperparameter of the gamma distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @seealso \code{\link{samplePower}}, \code{\link{samplePowerEst}}, \code{\link{find_size}}
#' @export samplePowerPoisson
#' @examples
#' \dontrun{
#' 
#' a1 <- 2; b1 <- 6
#' a2 <- 3; b2 <- 3
#' plot_gamma(c(a1, a2), c(b1, b2))
#' samplePowerPoisson(t = 10, a1, b1, a2, b2)
#' samplePowerPoisson(t = 11, a1, b1, a2, b2)
#' samplePowerPoisson(t = 10:11, a1, b1, a2, b2)
#' samplePowerPoisson(t = 250, a1, b1, a2, b2) # 56s
#' samplePowerEst(t = c(250, 500, 1000), a1, b1, a2, b2) # 5s, off by .001
#' 
#' a1 <- 4; b1 <- 2
#' a2 <- 6; b2 <- 2
#' plot_gamma(c(a1, a2), c(b1, b2))
#' samplePowerPoisson(t = 10, a1, b1, a2, b2)
#' 
#' 
#' }
samplePowerPoisson <- function(t, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1)
{

  # vectorize for multiple t
  if(is.numeric(t) && length(t) > 1){
    return(sapply(as.list(t), function(x){
      message(".", appendLF = FALSE)
      samplePowerPoisson(x, a1, b1, a2, b2, a, b, pi0, pi1, c)
    }))
  }

  # Rcpp
  samplePowerPoissonCpp(t, a1, b1, a2, b2, a, b, pi0, pi1, c)
}