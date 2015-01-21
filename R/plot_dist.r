#' Plot a beta density
#'
#' Plot a beta density
#' 
#' plot_beta creates a basic ggplot which can be styled.
#' 
#' @param alpha the beta's alpha parameter
#' @param beta the beta's beta parameter
#' @param n number of points to plot
#' @param ... ...
#' @return a ggplot object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @export plot_beta
#' @examples
#' 
#' \dontrun{
#'
#' plot_beta(1, 1)
#' plot_beta(.5, .5)
#' plot_beta(2, 2)
#' plot_beta(2, 2:3)
#'
#' # nonuniqueness of modes
#' plot_beta(4.5, 7.5) + ylim(0, 3)
#' plot_beta(1.12, 1.22) + ylim(0, 3)
#' plot_beta(c(4.5, 1.12), c(7.5, 1.22)) + ylim(0, 3)
#' 
#' 
#' plot_beta(.5, .5) + theme_bw()
#' plot_beta(.5, .5) + theme_classic()
#' 
#' theme_set(theme_bw())
#' plot_beta(.5, .5) + labs(x = "theta", y = "")
#' plot_beta(.5, .5) + ggtitle("my plot")
#' plot_beta(.5, .5) + labs(x = expression(theta[1]), y = "belief")
#' plot_beta(.5, .5) + labs(x = expression(paste(theta[1], "= 5")), y = "belief")
#' plot_beta(.5, .5) + labs(x = expression(paste(theta[1], "= ", 5)), y = "belief")
#'
#' }
#' 
plot_beta <- function(alpha, beta, n = 251, ...){
	
  # fake R CMD check
  param <- NULL; rm(param);	

  stopifnot(is.numeric(alpha))
  stopifnot(is.numeric(beta))  

  x <- seq(0, 1, length.out = n)
  x <- c(-.0001, x, 1.0001)

  if(length(alpha) == 1 && length(beta) == 1){  	
  	
    fx <- dbeta(x, alpha, beta)
    df <- data.frame(x = x, fx = fx)
  
    p <- ggplot(aes(x = x, y = fx), data = df) +
      geom_polygon(colour = "gray20", alpha = .25) +
      labs(y = 'likelihood')
      
  } else {
  	
    params <- as.data.frame(cbind(a = alpha, b = beta))
    fx <- as.data.frame( apply(params, 1, function(v) dbeta(x, v[1], v[2])) ) 
    combos <- apply(params, 1, function(v) paste0("alpha = ", v[1], " beta = ", v[2]))    
    names(fx) <- combos

    df <- cbind(x = x, fx)
    mdf <- melt(df, id = "x", var = "param", value.name = "fx")
    mdf$param <- as.character(mdf$param)    

    s <- apply(params, 1, function(v) 
      substitute(
        paste(alpha, " = ", x, ", ", beta, " = ", y),
        list(x = v[1], y = v[2]) 
      )
    )

    p <- ggplot(aes(x = x, y = fx, colour = param, fill = param), data = mdf) +
      geom_polygon(alpha = .25) +
      scale_colour_discrete(guide = "none") +
      scale_fill_discrete(labels = s) +
      labs(y = 'likelihood', colour = "", fill = "") +
      theme(legend.position = "top")  
  }
      
  p
}





















#' Plot a gamma density
#'
#' Plot a gamma density
#' 
#' plot_gamma creates a basic ggplot which can be styled.
#' 
#' @param alpha the gamma's alpha parameter
#' @param beta the gamma's beta parameter
#' @param n number of points to plot
#' @param q quantile used for xmax
#' @param ... ...
#' @return a ggplot object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @export plot_gamma
#' @examples
#' 
#' \dontrun{
#'
#' plot_gamma(1, 1)
#' plot_gamma(.5, .5)
#' plot_gamma(2, 2)
#' plot_gamma(1:5, 1)
#' plot_gamma(1, 1:5)
#'
#' plot_gamma(.5, .5, q = .999)
#' plot_gamma(.5, .5, q = .99)
#' plot_gamma(.5, .5, q = .9)
#'
# 
#' 
#' plot_gamma(2, 2) + theme_bw()
#' plot_gamma(2, 2) + theme_classic()
#' plot_gamma(2, c(2, 4)) + theme_classic()
#' 
#' theme_set(theme_bw())
#' plot_gamma(2, 2) + labs(x = "theta", y = "")
#' plot_gamma(2, 2) + ggtitle("my plot")
#' plot_gamma(2, 2) + labs(x = expression(lambda[1]), y = "belief")
#'
#' }
#' 
plot_gamma <- function(alpha, beta, n = 251, q = .999, ...){

  # fake R CMD check
  param <- NULL; rm(param);


  stopifnot(is.numeric(alpha))
  stopifnot(is.numeric(beta))  
   

  if(length(alpha) == 1 && length(beta) == 1){  	  	
  	
    xmax <- 1.05 * qgamma(q, shape = alpha, scale = beta)
  
    x <- seq(1E-5, xmax, length.out = n)
    x <- c(x[1], x, x[n])
    fx <- dgamma(x, shape = alpha, scale = beta)
    fx[1] <- fx[n+2] <- 0
    df <- data.frame(x = x, fx = fx)
  
    p <- ggplot(aes(x = x, y = fx), data = df) +
      geom_polygon(colour = "gray20", alpha = .25) +
      labs(y = 'likelihood')
      
  } else {
  	
    params <- as.data.frame(cbind(a = alpha, b = beta))
    
    xmax <- max(apply(params, 1, function(v){
      1.05 * qgamma(q, shape = v[1], scale = v[2])
    })) 
    
    x <- seq(1E-5, xmax, length.out = n)
    x <- c(x[1], x, x[n])    
    
    
    fx <- as.data.frame(apply(params, 1, function(v){
      y <- dgamma(x, shape = v[1], scale = v[2])
      y[1] <- y[n+2] <- 0
      y
    })) 
    combos <- apply(params, 1, function(v) paste0("alpha = ", v[1], " beta = ", v[2]))    
    names(fx) <- combos

    df <- cbind(x = x, fx)
    mdf <- melt(df, id = "x", var = "param", value.name = "fx")
    mdf$param <- as.character(mdf$param)    

    s <- apply(params, 1, function(v) 
      substitute(
        paste(alpha, " = ", x, ", ", beta, " = ", y),
        list(x = unname(v[1]), y = unname(v[2]))
      )
    )

    p <- ggplot(aes(x = x, y = fx, colour = param, fill = param), data = mdf) +
      geom_polygon(alpha = .25) +
      scale_colour_discrete(guide = "none") +
      scale_fill_discrete(labels = s) +
      labs(y = 'likelihood', colour = "", fill = "") +
      theme(legend.position = "top") 

  }
      
  p
}