#' Determine the probability of type I error of a test with given sample size
#'
#' Determine the probability of type I error of a test with given sample size
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
#' @seealso \code{\link{sampleAlphaBinomial}}, \code{\link{sampleAlphaPoisson}}, \code{\link{findSize}}
#' @export sampleAlpha
#' @examples
#' \dontrun{
#' 
#' 
#'  
#' 
#' # generate samples of size 30
#' N <- 1000000
#' n <- 30
#' 
#' a1 <- 3; b1 <- 7
#' a2 <- 7; b2 <- 3
#' plot_beta(c(a1, a2), c(b1, b2))
#' 
#' pi <- .3
#' pi <- rbeta(N, a1, b1)
#' data <- data.frame(
#'   x1 = rbinom(N, n, pi),
#'   x2 = rbinom(N, n, pi)  
#' )
#' test <- function(x) bayes_binom_test(x, n, a1, b1, a2, b2)$reject
#' mean( apply(data, 1, test) ) # .061832 at 1E6
#' sampleAlpha_ryan("binomial", n, 1, a1, b1) # .0204
#' sampleAlpha(n, a1, b1, a2, b2) 
#' 
#' 
#' 
#' 
#' # note that the significance depends on all parameters
#' a1 <- 3; b1 <- 7
#' a2 <- 90; b2 <- 10
#' plot_beta(c(a1, a2), c(b1, b2))
#' f <- function(x) bayes_binom_test(x, n, a1, b1, a2, b2)$reject
#' f(c(11, 8))
#' mean( apply(data, 1, f) )
#' sampleAlpha(n, a1, b1, a2, b2) 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' if(FALSE){
#' 
#' 
#' sim_alpha <- function(pi, N = 5E4, n = 30){
#' 
#'   data <- data.frame(
#'     x1 = rbinom(N, n, pi),
#'     x2 = rbinom(N, n, pi)  
#'   )
#' 
#' 
#'   a1 <- 3; b1 <- 7
#'   a2 <- 7; b2 <- 3
#' 
#'   f <- function(x) bayes_binom_test(x, n, a1, b1, a2, b2)$reject
#' 
#'   mean( apply(data, 1, f) )
#' }
#' 
#' sim_alpha(.3)
#' 
#' s <- seq(.01, .99, length.out = 250)
#' t1 <- sapply(as.list(s), function(x){
#'   message(".", appendLF = F)
#'   sim_alpha(x)
#' })
#' 
#' qplot(s, t1, geom = "line") +
#'   labs(x = expression(paste(pi, " = ", pi[1], " = ", pi[2])), y = expression(alpha)) +
#'   coord_equal()
#'   
#'   
#' sim_alpha_classic <- function(pi, N = 5E3, n = 30){
#' 
#'   data <- data.frame(
#'     x1 = rbinom(N, n, pi),
#'     x2 = rbinom(N, n, pi)  
#'   )
#'   
#'   f <- function(x) prop.test(x, c(n, n))$p.value < .05
#' 
#'   mean( apply(data, 1, f) )
#' }
#' sim_alpha_classic(.3)  
#'   
#' s <- seq(.01, .99, length.out = 25)
#' t1 <- sapply(as.list(s), function(x){
#'   message(".", appendLF = F)
#'   sim_alpha_classic(x)
#' })
#' 
#' qplot(s, t1, geom = "line") +
#'   labs(x = expression(paste(pi, " = ", pi[1], " = ", pi[2])), y = expression(alpha))
#' 
#' 
#' 
#' 
#' 
#' }
#' 
#' 
#'
#' }
sampleAlpha <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1,
  family = c("binomial", "poisson"))
{

  # check arguments
  stop2 <- function(x) stop(x, call. = FALSE)
  if(missing(n)) stop2("the sample size n must be specified.  see ?sample_power.")
  if(missing(a1)) stop2("the hyperparameter a1 must be specified.  see ?sample_power.")  
  if(missing(b1)) stop2("the hyperparameter b1 must be specified.  see ?sample_power.")  
  if(missing(a2)) stop2("the hyperparameter a2 must be specified.  see ?sample_power.")  
  if(missing(b2)) stop2("the hyperparameter b2 must be specified.  see ?sample_power.")        
  family <- match.arg(family)
  stopifnot( pi0 >= 0 )
  stopifnot( pi1 >= 0 )  
  stopifnot( pi0 + pi1 == 1 )

  # dispatch correct function
  if(family == "binomial"){
    out <- sampleAlphaBinomial(n, a1, b1, a2, b2, a, b, pi0, pi1, c)
  }

  if(family == "poisson"){
    out <- sampleAlphaPoisson(n, a1, b1, a2, b2, a, b, pi0, pi1, c)
  }

  # return
  out
}


















































#' Compute the probability of type I error of the binomial proportions test for a given sample size (n) with specified beta priors
#'
#' Compute the probability of type I error of the binomial proportions test for a given sample size (n) with specified beta priors
#' 
#' @param n the sample size 
#' @param a1 alpha, the hyperparameter of the prior distribution
#' @param b1 beta, the hyperparameter of the prior distribution
#' @param a2 alpha, the hyperparameter of the prior distribution
#' @param b2 beta, the hyperparameter of the prior distribution
#' @param a alpha, the hyperparameter of the beta distribution under the null
#' @param b beta, the hyperparameter of the beta distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @seealso \code{\link{sampleAlpha}}, \code{\link{findSize}}
#' @export sampleAlphaBinomial
#' @examples
#' \dontrun{
#' 
#' a <- 5; b <- 5
#' sampleAlpha(t = 20, a, b)
#' 
#' 
#' }
sampleAlphaBinomial <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1)
{

  # vectorize for multiple n
  if(is.numeric(n) && length(n) > 1){
    return(sapply(as.list(n), function(x){
      message(".", appendLF = FALSE)
      sampleAlphaBinomial(x, a1, b1, a2, b2, a, b, pi0, pi1, c)
    }))
  }
  
  # enumerate possible outcomes - this gets huge
  df <- expand.grid(y1 = 0:n, y2 = 0:n)
  
  # test, TRUE when fail to reject
  test <- function(y1, y2){
    lbeta(y1 + y2 + a, 2 * n - y1 - y2 + b) - lbeta(a, b) -
      lbeta(y1 + a1, n - y1 + b1) + lbeta(a1, b1) -
      lbeta(y2 + a2, n - y2 + b2) + lbeta(a2, b2) > log(c * pi1/pi0)
  }
  
  # and put its probability here between the lines return(exp(... and ...)) 
  value <- function(y1, y2){
  	if(test(y1, y2)){  	
      return(
        exp(lchoose(n, y1) + lbeta(y1 + a1, n - y1 + b1) -
          lbeta(a1, b1) + lchoose(n, y2) + lbeta(y2 + a2,
          n - y2 + b2) - lbeta(a2, b2))
      )
    } else {
      return(0)
    }
  }
  
  # sum the probabilities of each of the rejected y's and return
  1 - sum(apply(df, 1, function(v) value(v[1], v[2])))
}










































#' Compute the probability of type I error of the Poisson rate test for a given sample size (t) with specified gamma priors
#'
#' Compute the probability of type I error of the Poisson rate test for a given sample size (t) with specified gamma priors
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
#' @seealso \code{\link{sampleAlpha}}, \code{\link{findSize}}
#' @export sampleAlphaPoisson
#' @examples
#' \dontrun{
#' 
#' a <- ; b <- 
#' sampleAlpha(t = 20, a, b)
#' 
#' 
#' }
sampleAlphaPoisson <- function(t, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1)
{

  # vectorize for multiple t
  if(is.numeric(t) && length(t) > 1){
    return(sapply(as.list(t), function(x){
      message(".", appendLF = FALSE)
      sampleAlphaPoisson(x, a1, b1, a2, b2, a, b, pi0, pi1, c)
    }))
  }

  # determine essential support of the poisson (/gamma mixture)
  y1min <- qnbinom(.00001, size = a1, prob = b1/(t+b1))
  y2min <- qnbinom(.00001, size = a2, prob = b2/(t+b2))   
	
  y1max <- qnbinom(.99999, size = a1, prob = b1/(t+b1))
  y2max <- qnbinom(.99999, size = a2, prob = b2/(t+b2))  
  
  # enumerate possible outcomes - this gets huge
  df <- expand.grid(y1 = y1min:y1max, y2 = y2min:y2max)
  
  # test, TRUE when fail to reject
  test <- function(y1, y2){
    a*log(b) + lgamma(y1+y2+a) - lgamma(a) - (y1+y2+a)*log(2*t+b) + 
      lgamma(a1) + lgamma(a2) + (y1+a1)*log(t+b1) + (y2+a2)*log(t+b2) - 
      a1*log(b1) - a2*log(b2) - lgamma(y1+a1) - lgamma(y2+a2) > log(c * pi1/pi0)    
  }
  
  # and put its probability here between the lines return(exp(... and ...)) 
  value <- function(y1, y2){
  	if(test(y1, y2)){  	
      return(
        exp( 
          (y1+y2)*log(t) + a1*log(b1) + a2*log(b2) + lgamma(y1+a1) +
          lgamma(y2+a2) - lfactorial(y1) - lfactorial(y2) - lgamma(a1) - lgamma(a2) -
          (y1+a1)*log(t+b1) - (y2+a2)*log(t+b2) 
        )
      )
    } else {
      return(0)
    }
  }
  
  # sum the probabilities of each of the rejected y's and return
  1 - sum(apply(df, 1, function(v) value(v[1], v[2])))
}