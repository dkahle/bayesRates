#' Make a power curve
#'
#' Make a power curve
#' 
#' @param ts vector specifying x-range
#' @param a1 alpha, the hyperparameter of the gamma distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the gamma distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the gamma distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the gamma distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the gamma distribution under the null
#' @param b beta, the hyperparameter of the gamma distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param family "binomial" or "poisson", depending on test
#' @param ylim y limits
#' @export plotPower
#' @examples
#' \dontrun{
#' 
#' a1 <- 2; b1 <- 6
#' a2 <- 3; b2 <- 3
#' plotPower(ts = c(0,250), a1, b1, a2, b2)
#' plotPower(ts = c(0,1000), a1, b1, a2, b2)
#' plotPower(ts = c(0,100), a1, b1, a2, b2)
#' plotPower(ts = 0:300, a1, b1, a2, b2, method = 'exact')
#' # this last plot shows that the approximation is not quite exact for small values
#' # but decent for larger ones
#'
#' # styling, all ggplot2 styling works
#' plotPower(t = c(0,1000), a1, b1, a2, b2) +
#'   theme_bw() + opts(title = 'My Power Plot')
#' 
#' 
#' 
#' a1 <- 4; b1 <- 4
#' a2 <- 8; b2 <- 4
#' plotPower(t = c(0,200), a1, b1, a2, b2, family = "poisson") +
#'   geom_hline(yintercept = .80)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' a1 <- 4; b1 <- 2
#' a2 <- 6; b2 <- 2
#' plotPower(t = c(0,1000), a1, b1, a2, b2)
#' plotPower(t = c(0,25), a1, b1, a2, b2)
#' plotPower(t = 0:25, a1, b1, a2, b2, method = 'exact')
#' 
#' 
#' 
#' 
#' a1 <- c(2,4); b1 <- c(6,2);
#' a2 <- c(3,6); b2 <- c(3,2);
#' plotPower(t = c(0,1000), a1, b1, a2, b2)
#'
#' library(ggplot2)
#' last_plot() + theme_bw()
#' 
#' 
#' 
#' 
#'  }
plotPower <- function(ts, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c = 1,
  family = c("binomial", "poisson"),  ylim = 0:1)
{   
  est_power <- NULL; rm(est_power);
  
  # check args
  family <- match.arg(family)
  if(min(ts) == 0) ts[ts == 0] <- 1L
  	
  # is a range given (length 2), or a set of t's?
  if(length(ts) == 2){
    t <- seq(min(ts), max(ts), 1)
  } else {
    t <- ts
  }
  	
  # t range is already provided
  df <- data.frame(
    t = t,
    power = samplePower(t, a1, b1, a2, b2, a, b, pi0, pi1, c, family)
  )
    
    
  # remove 0's
  zeros <- with(df, which(power == 0))
  if(length(zeros) > 0) df <- df[-zeros,]
  nans <- with(df, which(is.nan(power)))
  if(length(nans) > 0) df <- df[-nans,]
  
    
  # model
  power_mod <- lm(power ~ log(t+1e-3), data = df)
  ts <- pmax(extendrange(range(ts)),0)
  ts <- seq(ts[1], ts[2], length.out = 501)
  df_line <- data.frame(
    t = ts, 
    power = pmin(predict(power_mod, data.frame(t = ts)),1)
  )    
    
    
  p <- ggplot(aes(x = t, y = power), data = df_line) +
    #geom_path() + 
    ylim(ylim) +
    geom_point(alpha = 2/4, colour = "red", data = df) +
    geom_line(alpha = 2/4, colour = "red", data = df) +
    #geom_hline(yintercept = 0, colour = "black", alpha = .5) +
    #geom_hline(yintercept = 1, colour = "black", alpha = .5) +    
    labs(x = "t (Sample Size)", y = "Power")         	  	
  
  p     
}