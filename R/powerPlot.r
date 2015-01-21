#' Make a power curve
#'
#' Make a power curve
#' 
#' @param t_range vector specifying x-range
#' @param a1 alpha, the hyperparameter of the gamma distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the gamma distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the gamma distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the gamma distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the gamma distribution under the null
#' @param b beta, the hyperparameter of the gamma distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param ylim y limits
#' @param method method argument to pass to \code{\link{samplePowerEst}}
#' @seealso \code{\link{samplePowerEst}}
#' @export powerPlot
#' @examples
#' \dontrun{
#' 
#' a1 <- 2; b1 <- 6
#' a2 <- 3; b2 <- 3
#' powerPlot(t = c(0,1000), a1, b1, a2, b2)
#' powerPlot(t = c(0,100), a1, b1, a2, b2)
#' powerPlot(t = 0:100, a1, b1, a2, b2, method = 'exact')
#' # this last plot shows that the approximation is not quite exact for small values
#' # but decent for larger ones
#'
#' # styling, all ggplot2 styling works
#' powerPlot(t = c(0,1000), a1, b1, a2, b2) +
#'   theme_bw() + opts(title = 'My Power Plot')
#' 
#' 
#' 
#' 
#' a1 <- 4; b1 <- 2
#' a2 <- 6; b2 <- 2
#' powerPlot(t = c(0,1000), a1, b1, a2, b2)
#' powerPlot(t = c(0,25), a1, b1, a2, b2)
#' powerPlot(t = 0:25, a1, b1, a2, b2, method = 'exact')
#' 
#' 
#' 
#' 
#' a1 <- c(2,4); b1 <- c(6,2);
#' a2 <- c(3,6); b2 <- c(3,2);
#' powerPlot(t = c(0,1000), a1, b1, a2, b2)
#'
#' library(ggplot2)
#' last_plot() + theme_bw()
#' 
#' 
#' 
#' 
#'  }
powerPlot <- function(t_range, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, ylim = 0:1, method = 'fastest')
{   
  est_power <- NULL; rm(est_power);
  
  # check args
  match.arg(method, c('fastest','fast','moderate','exact'))
  
  if(method != 'exact'){
  
    # make range of t's to approx for curve (note : not integer)
    xlim <- pmax(extendrange(t_range), 0)  
    x <- seq(xlim[1], xlim[2], length.out = 501)
  
    # approximate with samplePowerEst (log)
    df <- NULL
  
    for(k in seq_along(a1)){
      id <- paste0('a1 = ', a1[k], ', b1 = ', b1[k], '\n', 
  	    'a2 = ', a2[k], ', b2 = ', b2[k], '\n', 
  	    'a = ', a[k], ', b = ', b[k], '\n'
      )
  	
      df_tmp <- data.frame(
        samplePowerEst(x, a1 = a1[k], b1 = b1[k], a2 = a2[k], 
          b2 = b2[k], a = a[k], b = b[k], pi0 = pi0, pi1 = pi1, 
          method = method),
        id = id
      )
  
      # set 0 sample size 
      df_tmp$t[1] <- 0
      df_tmp$est_power[1] <- 0  
      df_tmp$time[1] <- 0        


      df <- rbind(df, df_tmp)
    }    
    
    # make / return plot  
    if(length(a1) == 1){
      p <- qplot(t, est_power, data = df, geom = 'path', ylim = ylim)
    } else {
      p <- qplot(t, est_power, data = df, colour = id, geom = 'path', ylim = ylim)
      p <- p + labs(colour = 'Parameter\nConfigurations')    
    }
  
    p <- p +
      geom_hline(yintercept = 0, colour = 'black', alpha = .5) +
      geom_hline(yintercept = 1, colour = 'black', alpha = .5) +    
      labs(x = 't (Sample Size)', y = 'Estimated Power') 
      
  } else if(method == 'exact'){
  	
  	# is a range given (length 2), or a set of t's?
    if(length(t_range) == 2){
      t <- seq(min(t_range), max(t_range), length.out = 8)
    } else {
      t <- t_range
    }
  	
    # t range is already provided
    df <- data.frame(
      t = t,
      power = samplePower(t, a1, b1, a2, b2, a, b, pi0, pi1)
    )
    
    
    # remove 0's
    zeros <- with(df, which(power == 0))
    if(length(zeros) > 0) df <- df[-zeros,]
    nans <- with(df, which(is.nan(power)))
    if(length(nans) > 0) df <- df[-nans,]
    
    
    # model
    power_mod <- lm(power ~ log(t), data = df)
    t_range <- pmax(extendrange(range(t_range)),0)
    ts <- seq(t_range[1], t_range[2], length.out = 501)
    df_line <- data.frame(
      t = ts, 
      power = pmin(predict(power_mod, data.frame(t = ts)),1)
    )    
    
    
    p <- ggplot(aes(x = t, y = power), data = df_line) +
      geom_path() + ylim(ylim) +
      geom_point(alpha = 2/4, colour = 'red', data = df) +
      geom_hline(yintercept = 0, colour = 'black', alpha = .5) +
      geom_hline(yintercept = 1, colour = 'black', alpha = .5) +    
      labs(x = 't (Sample Size)', y = 'Power')         	  	
  } 
  
  p     
}