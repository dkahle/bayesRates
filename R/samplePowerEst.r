#' Approximate the power of a given sample size
#'
#' Approximate a sample size
#' 
#' @param n the sample size - possibly a vector
#' @param a1 alpha, the hyperparameter of the gamma distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the gamma distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the gamma distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the gamma distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the gamma distribution under the null
#' @param b beta, the hyperparameter of the gamma distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param method "fastest", "fast", "moderate", decreasing in speed and increasing in accuracy
#' @param family "binomial" or "poisson"
#' @seealso \code{\link{samplePower}}, \code{\link{findSize}}
#' @export
#' @examples
#' \dontrun{
#' 
#' # binomial
#' a1 <- 2; b1 <- 6
#' a2 <- 3; b2 <- 3
#' samplePower(n = 10, a1, b1, a2, b2)
#' samplePowerEst(n = 10, a1, b1, a2, b2)
#' 
#' # poisson
#' a1 <- 4; b1 <- 2
#' a2 <- 6; b2 <- 2
#' samplePower(n= 20, a1, b1, a2, b2, family = "poisson")
#' samplePowerEst(n = 20, a1, b1, a2, b2, family = "poisson")
#' samplePowerEst(n = 20, a1, b1, a2, b2, family = "poisson", method = "fast")
#' samplePowerEst(n = 20, a1, b1, a2, b2, family = "poisson", method = "moderate")
#' 
#' # this takes a moment
#' samplePowerEst(n = 100*1:10, a1, b1, a2, b2, family = "poisson")
#' samplePower(n = 200, a1, b1, a2, b2, family = "poisson")
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
#' ## predict power for a given sample size (poisson)
#' ############################################################
#' n <- 100
#' mult <- 1
#' power <- numeric(n)
#' time <- numeric(n)
#' t <- c(mult*(1:n),150,200,250,300,350,400,450,500,750,1000)
#' t <- c(1, 50, 100 ,150,200,250,300,350,400,450,500,750,1000)
#' for(k in seq_along(t)){
#'   start <- Sys.time()
#'   message(paste0("t = ", t[k], "... "), appendLF = FALSE)   
#'   power[k] <- samplePower(t[k], a1, b1, a2, b2, family = "poisson")
#'   time[k] <- end <- difftime(Sys.time(), start, units = "secs")
#'   if(k == length(t)) message("")
#' }
#' df <- data.frame(n = t, time = time, power = power)
#' 
#' 
#' library(ggplot2)
#' # investigate power as a function of sample size n
#' qplot(n, power, data = df) + ylim(0,1)
#' 
#' qplot(n, power, data = df) +
#'   stat_smooth(method = "lm", formula = y ~ log(x)) +
#'   ylim(0,1)
#'   
#' mod <- lm(power ~ log(n), data = df)  
#' s <- 1:1000
#' pdf <- data.frame(
#'   n = s,
#'   power = predict(mod, data.frame(n = s))
#' )  
#' qplot(n, power, data = pdf, geom = "path", ylim = c(0,1)) +
#'     labs(x = "n (Sample Size)", y = "Estimated Power") +
#'   geom_point(aes(x = n, y = power), data = df, color = "red")
#' 
#' 
#' 
#' ## predict time necessary for the exact power to be computed as a 
#' ## function of the sample size (poisson)
#' ############################################################
#' 
#' qplot(n, time, data = df)
#' 
#' # not log...
#' qplot(n, time, data = df) + scale_y_log10()
#' 
#' # quadratic
#' qplot(n, time, data = df) +
#'   stat_smooth(method = "lm", formula = y ~ x + I(x^2))
#'   
#' mod <- lm(time ~ n + I(n^2), data = df)
#' s <- 1:1000
#' tdf <- data.frame(
#'   n = s,
#'   time = predict(mod, data.frame(n = s))
#' )  
#' qplot(n, time, data = tdf, geom = "path") +
#'   labs(x = "n (Sample Size)", y = "Time (Seconds)")
#' 
#' 
#' start <- Sys.time()
#' samplePower(500, a1, b1, a2, b2, family = "poisson")
#' difftime(Sys.time(), start, units = "secs")
#' subset(tdf, n == 500)
#' # predicted times are very close
#' 
#' 
#' 
#' }
samplePowerEst <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, 
  method = c("fastest", "fast", "moderate"), 
  family = c("binomial", "poisson")
){
  
  # check arguments
  method <- match.arg(method)
  family <- match.arg(family)

  if(method == "fastest"){ 
    n <- 10
    mult <- 4 
  } else if(method == "fast") {
    n <- 20
    mult <- 8
  } else { # moderate
    n <- 40
    mult <- 30
  }
  
  sample_points <- mult * 1:n
  
  # redistribute sample points if they exceed max(n)
  if(max(n) < max(sample_points)){
  	sample_points <- unique(round(seq(1, max(n), length.out = n)))
    n <- length(sample_points)
  }
  
  # calculate exact power of sample points
  power <- numeric(n)
  time <- numeric(n)
  for(k in 1:n){
    start <- Sys.time()
    power[k] <- samplePower(sample_points[k], a1, b1, a2, b2, family = family)
    time[k] <- end <- difftime(Sys.time(), start, units = "secs")
    message("..", appendLF = FALSE)
  }
  cat("\n")

  
  # arrange df
  df <- data.frame(n = sample_points, time = time, power = power)  
  
  # make models
  power_mod <- lm(power ~ log(n), data = df)  
  time_mod <- lm(time ~ n + I(n^2), data = df)    

  # predict and return data.frame
  data.frame(
    n = n,
    est_power = pmin(predict(power_mod, data.frame(n = n)), 1),
    est_time = predict(time_mod, data.frame(n = n))
  )    
}