#' Bayesian two-sample test of proportions
#'
#' Bayesian two-sample test of proportions
#' 
#' @param x  two-vector containing successful trials
#' @param n  number of trials (the same for both samples)
#' @param a1 alpha, the hyperparameter of the beta distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the beta distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the beta distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the beta distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the beta distribution under the null
#' @param b beta, the hyperparameter of the beta distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c1 loss associated with type I error
#' @param c2 loss associated with type I error
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param rule the decision rule for the Bayes factor B.  by default, test rejects when B < pi1/pi0 * c2/c1
#' @seealso \code{\link{sample_power}}, \code{\link{sample_power_est}}, \code{\link{find_size}}, \code{\link{plot_binom_rule}}
#' @export bayes_binom_test
#' @examples
#' \dontrun{
#' 
#' a1 <- 3/5; b1 <- 7/5
#' a2 <- 7/5; b2 <- 3/5
#' plot_beta(c(a1, a2), c(b1, b2))
#' bayes_binom_test(c(12, 16), 20, a1, b1, a2, b2) # reject
#' bayes_binom_test(c(13, 16), 20, a1, b1, a2, b2) # fail to reject
#'
#' # compare to :
#' prop.test(c(9, 16), c(20, 20)) # almost exact reject at 5% level
#' prop.test(c(10, 16), c(20, 20)) # fail to reject
#' prop.test(c(11, 16), c(20, 20)) # fail to reject
#'
#'
#' # but even slightly more informative priors can overpower results
#' # watch the proportions to see why
#' a1 <- 3; b1 <- 7
#' a2 <- 7; b2 <- 3
#' plot_beta(c(a1, a2), c(b1, b2))
#' bayes_binom_test(c(12, 16), 20, a1, b1, a2, b2) # reject
#' bayes_binom_test(c(13, 16), 20, a1, b1, a2, b2) # reject
#' bayes_binom_test(c(16, 16), 20, a1, b1, a2, b2) # reject (!)
#' bayes_binom_test(c(17, 16), 20, a1, b1, a2, b2) # reject (!)
#' bayes_binom_test(c(18, 16), 20, a1, b1, a2, b2) # fail to reject (!)
#' 
#' 
#' }
bayes_binom_test <- function(x, n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c1 = 1, c2 = 1, c = c1/c2, rule = pi1/pi0 * c2/c1)
{
  
  x1 <- x[1]; x2 <- x[2]
  
  logPosteriorOddsNum <- lchoose(n, x1) + lchoose(n, x2) +
    dbetabinom(x1 + x2, 2*n, a, b, log = TRUE) + log(pi0)  
  
  logPosteriorOddsDenom <- lchoose(2*n, x1 + x2) +
    dbetabinom(x1, n, a1, b1, log = TRUE) + dbetabinom(x2, n, a2, b2, log = TRUE) +
    log(pi1)
    
  logPosteriorOdds <- logPosteriorOddsNum - logPosteriorOddsDenom
  
  B <- exp( logPosteriorOdds - log(pi0) + log(pi1) )
  attr(B, "names") <- "Bayes factor"

  px <- capture.output(print(dput(unname(x))))[1]
  pn <- capture.output(print(dput(rep(n, 2))))[1]
  
  ests <- c("prop 1" = (a1 + x1) / (a1 + b1 + n), "prop 2" = (a2 + x2) / (a2 + b2 + n))

  params <- matrix(c(a, a1, a2, b, b1, b2), nrow = 3, ncol = 2,
    dimnames = list(
    c("null", "alt 1", "alt 2"), 
    c("alpha", "beta"))
  )
  
  if(B <= rule){
    conclusion <- "Conclusion : Reject null hypothesis of proportion equivalence"
    reject <- TRUE
  } else {
    conclusion <- "Conclusion : Fail to reject null hypothesis of proportion equivalence"
    reject <- FALSE
  }

  out <- list(
    method = "2-sample test for equality of proportions with prior information",
    statistic = B,
    odds = exp(logPosteriorOdds),
    estimate = ests,
    data.name = paste(px, "out of", pn),
    rule = c("Bayes rule critical value" = rule),
    params = params,
    conclusion = conclusion,
    hyp_priors = c("P(H0)" = pi0, "P(H1)" = pi1),
    reject = reject
  )
  class(out) <- "btest"

  out
}














#' Density of the beta-binomial distribution with parameters size, alpha, and beta
#'
#' Density of the beta-binomial distribution with parameters size, alpha, and beta
#' 
#' @param x observation
#' @param size the beta-binomial size parameter
#' @param alpha the beta-binomial alpha parameter
#' @param beta the beta-binomial beta parameter
#' @param log return log value?
#' @export dbetabinom
#' @examples
#' \dontrun{
#' 
#' dbetabinom(2, 10, .5, .5)
#' dbetabinom(0:10, 10, .5, .5)
#' sum(dbetabinom(0:10, 10, .5, .5))
#' 
#' }
dbetabinom <- function(x, size, alpha, beta, log = FALSE){
	
  if(length(x) > 1) return(
    sapply(as.list(x), function(.) dbetabinom(., size, alpha, beta, log = log))
  )
	
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  stopifnot(x >= 0)	
  stopifnot(is.wholenumber(x))    
  stopifnot(size > 0)	
  stopifnot(is.wholenumber(size))  
  stopifnot(alpha > 0)
  stopifnot(beta > 0)  

  logP <- lchoose(size, x) + lbeta(x + alpha, size - x + beta) - lbeta(alpha, beta)  
  if(log) return(logP)

  exp(logP)
}







































#' Bayesian two-sample test of Poisson rates
#'
#' Bayesian two-sample test of Poisson rates
#' 
#' @param x  two-vector containing occurrences
#' @param t  duration of trials (the same for both samples)
#' @param a1 alpha, the hyperparameter of the gamma distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the gamma distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the gamma distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the gamma distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the gamma distribution under the null
#' @param b beta, the hyperparameter of the gamma distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param rule the decision rule for the Bayes factor B.  by default, test rejects when B < pi1/pi0 * c2/c1
#' @seealso \code{\link{sample_power}}, \code{\link{sample_power_est}}, \code{\link{find_size}}
#' @export bayes_pois_test
#' @examples
#' \dontrun{
#' 
#' a1 <- 5; b1 <- .2 # E = 1
#' a2 <- 30; b2 <- .167 # E = 5
#' plot_gamma(c(a1, a2), c(b1, b2))
#' bayes_pois_test(c(18, 105), 20, a1, b1, a2, b2)
#' bayes_pois_test(c(18, 165), 20, a1, b1, a2, b2)
#'
#' # compare to :
#' prop.test(c(10, 16), c(20, 20))
#' 
#' 
#'  }
bayes_pois_test <- function(x, t, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1, rule = c*pi1/pi0)
{
  
  x1 <- x[1]; x2 <- x[2]

  logB <- a*log(b) + lgamma(x1+x2+a) - lgamma(a) - (x1+x2+a)*log(2*t+b) + 
    lgamma(a1) + lgamma(a2) + (x1+a1)*log(t+b1) + (x2+a2)*log(t+b2) - 
    a1*log(b1) - a2*log(b2) - lgamma(x1+a1) - lgamma(x2+a2)
  attr(logB, "names") <- "Bayes factor"

  px <- capture.output(print(dput(unname(x))))[1]
  pn <- capture.output(print(dput(rep(t, 2))))[1]
  
  ests <- c("rate 1" = (a1 + x1) / (b1 + t), "rate 2" = (a2 + x2) / (b2 + t))

  params <- matrix(c(a, a1, a2, b, b1, b2), nrow = 3, ncol = 2,
    dimnames = list(
    c("null", "alt 1", "alt 2"), 
    c("alpha", "beta"))
  )
  
  if(logB <= log(rule)){
    conclusion <- "Conclusion : Reject null hypothesis of rate equivalence"
    reject <- TRUE
  } else {
    conclusion <- "Conclusion : Fail to reject null hypothesis of rate equivalence"
    reject <- FALSE    
  }

  out <- list(
    method = "2-sample test for equality of Poisson rates with prior information",
    statistic = exp(logB),
    estimate = ests,
    data.name = paste(px, "out of", pn),
    rule = c("Bayes rule critical value" = rule),
    params = params,
    conclusion = conclusion,
    hyp_priors = c("P(H0)" = pi0, "P(H1)" = pi1),
    reject = reject
  )
  class(out) <- "btest"

  out
}










































#' Print a Bayesian test
#'
#' Print a Bayesian test
#' 
#' @param x an object of class btest
#' @param digits number of digits to show
#' @param ... additional parameters
#' @usage \method{print}{btest}(x, digits = 4L, ...)
#' @return invisible object
#' @export
#' @examples
#' \dontrun{
#'   bayes_binom_test(11:12, 20, 1, 1, 1, 1)
#' }
#' 
print.btest <- function (x, digits = 4L, ...){

    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")

    out <- character()
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 
            4))))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")

    if (!is.null(x$hyp_priors)) {
        cat("Hypotheses priors:\n")
        print(x$hyp_priors, ...)
    }    
    
    if (!is.null(x$hyp_priors)) {
        cat("Prior odds:\n")
        podds <- unname(x$hyp_priors)
        print(podds[1] / podds[2], ...)
    }        

    if (!is.null(x$params)) {
        cat("Prior hyper-parameters:\n")
        print(x$params, ...)
    }    

    if (!is.null(x$estimate)) {
        cat("Sample estimates:\n")
        print(x$estimate, ...)
    }
    
    if (!is.null(x$odds)) {
        cat("Posterior odds:\n")
        print(x$odds, ...)
    }        
    
    out <- character()
    if (!is.null(x$rule)) 
        out <- c(out, paste(names(x$rule), "=", format(round(x$rule, 
            4))))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    
    cat(strwrap(x$conclusion), sep = "\n")    
    
    cat("\n")
    
    invisible(x)
}