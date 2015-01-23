#' Plot a decision rule
#'
#' Plot a decision rule
#' 
#' @param n  number of trials (the same for both samples)
#' @param a1 alpha, the hyperparameter of the prior distribution of the first poisson rate
#' @param b1 beta, the hyperparameter of the prior distribution of the first poisson rate
#' @param a2 alpha, the hyperparameter of the prior distribution of the second poisson rate
#' @param b2 beta, the hyperparameter of the prior distribution of the second poisson rate
#' @param a alpha, the hyperparameter of the prior distribution under the null
#' @param b beta, the hyperparameter of the prior distribution under the null
#' @param pi0 the prior probability of the null hypothesis
#' @param pi1 the prior probability of the alternative hypothesis
#' @param c1 loss associated with type I error
#' @param c2 loss associated with type I error
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param rule the decision rule for the Bayes factor B.  by default, test rejects when B < pi1/pi0 * c2/c1
#' @param sizeby size the outcomes by the probability under the null or alternative
#' @param geom ggplot2 geom used
#' @return a ggplot object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @export plotRule
#' @examples
#' 
#' \dontrun{
#'  
#' plotBinomRule(30, 1, 1, 1, 1)
#' plotBinomRule(30, 3, 7, 7, 3, sizeby = "null") + theme_bw()
#' plotBinomRule(30, 3, 7, 7, 3, sizeby = "alt") + theme_bw()
#'
#' }
#' 
plotRule <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c1 = 1, c2 = 1, c = c1/c2, rule = pi1/pi0 * c2/c1,
  sizeby = c("none", "null", "alt"), geom = c("auto", "point", "tile"))
{
  
  # fake R CMD check
  y1 <- NULL; rm(y1);
  y2 <- NULL; rm(y2);  
  p <- NULL; rm(p);    
  action <- NULL; rm(action);  
  sizeby <- match.arg(sizeby)
  geom <- match.arg(geom)  
  
  if(geom == "auto" && n <= 50){
    geom <- "point"
  } else if(geom == "auto" && n > 50){
    geom <- "tile"  	
  }
  
  if(sizeby != "none") geom <- "point"
  
  df <- expand.grid(y1 = 0:n, y2 = 0:n)
  df$action <- apply(df, 1, function(x) 
    bayesBinomTest(x, n, a1, b1, a2, b2, a, b, pi0, pi1, c1, c2, c, rule)$reject
  )
  
  df$action[df$action == TRUE] <- "Reject"
  df$action[df$action == "FALSE"] <- "Accept"  
  
  if(geom == "tile") return(
    qplot(y1, y2, data = df, geom = "tile", fill = action) +
      scale_fill_manual("Decision",
        values = c(
          "Accept" = rgb(61, 153, 86, maxColorValue = 256),
          "Reject" = rgb(255, 10, 27, maxColorValue = 256)
        )
      )
  )
  
  if(sizeby == "none") return(    
    qplot(y1, y2, data = df, geom = "point", color = action, shape = I(15)) +
      scale_color_manual("Decision",
        values = c(
          "Accept" = rgb(61, 153, 86, maxColorValue = 256),
          "Reject" = rgb(255, 10, 27, maxColorValue = 256)
        )
      )    
  )
  
  if(sizeby == "null"){
    prob <- function(v) exp(
      lchoose(n, v[1]) + lchoose(n, v[2]) + dbetabinom(sum(v), 2*n, a, b, log = T) -
        lchoose(2*n, sum(v))
    )
  } else { # sizeby == "alt"
    prob <- function(v) exp(
      dbetabinom(v[1], n, a1, b1, log = T) + dbetabinom(v[2], n, a2, b2, log = T)
    )
  }
  
  df$p <- apply(df[,c("y1","y2")], 1, prob)
  
  qplot(y1, y2, data = df, geom = "point", color = action, shape = I(15), size = p) +
    scale_color_manual("Decision",
      values = c(
        "Accept" = rgb(61, 153, 86, maxColorValue = 256),
        "Reject" = rgb(255, 10, 27, maxColorValue = 256)
      )
    ) +
    scale_size("Probability")    
}






















#' Plot a binomial decision rule
#'
#' Plot a binomial decision rule
#' 
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
#' @param sizeby size the outcomes by the probability under the null or alternative
#' @param geom ggplot2 geom used
#' @return a ggplot object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @export plotBinomRule
#' @examples
#' 
#' \dontrun{
#'	
#' plotBinomRule(30, 1, 1, 1, 1)
#' plotBinomRule(100, 1, 1, 1, 1)
#' plotBinomRule(30, 1, 1, 1, 1, sizeby = "null")
#' plotBinomRule(30, 1, 1, 1, 1, sizeby = "alt") # uniform
#' 
#' plotBinomRule(30, 3, 7, 7, 3)
#' plotBinomRule(100, 3, 7, 7, 3)
#'
#' plotBinomRule(30, 3, 7, 7, 3, sizeby = "null") + theme_bw()
#' plotBinomRule(30, 3, 7, 7, 3, sizeby = "alt") + theme_bw()
#'
#' }
#' 
plotBinomRule <- function(n, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c1 = 1, c2 = 1, c = c1/c2, rule = pi1/pi0 * c2/c1,
  sizeby = c("none", "null", "alt"), geom = c("auto", "point", "tile"))
{

  # fake R CMD check
  y1 <- NULL; rm(y1);
  y2 <- NULL; rm(y2);  
  p <- NULL; rm(p);    
  action <- NULL; rm(action);  
  sizeby <- match.arg(sizeby)
  geom <- match.arg(geom)  
  
  if(geom == "auto" && n <= 50){
    geom <- "point"
  } else if(geom == "auto" && n > 50){
    geom <- "tile"  	
  }
  
  if(sizeby != "none") geom <- "point"
  
  df <- expand.grid(y1 = 0:n, y2 = 0:n)
  df$action <- apply(df, 1, function(x) 
    bayesBinomTest(x, n, a1, b1, a2, b2, a, b, pi0, pi1, c1, c2, c, rule)$reject
  )
  
  df$action[df$action == TRUE] <- "Reject"
  df$action[df$action == "FALSE"] <- "Accept"  
  
  if(geom == "tile") return(
    qplot(y1, y2, data = df, geom = "tile", fill = action) +
      scale_fill_manual("Decision",
        values = c(
          "Accept" = rgb(61, 153, 86, maxColorValue = 256),
          "Reject" = rgb(255, 10, 27, maxColorValue = 256)
        )
      )
    )
  
  if(sizeby == "none") return(    
    qplot(y1, y2, data = df, geom = "point", color = action, shape = I(15)) +
      scale_color_manual("Decision",
        values = c(
          "Accept" = rgb(61, 153, 86, maxColorValue = 256),
          "Reject" = rgb(255, 10, 27, maxColorValue = 256)
        )
      )    
    )
    
  if(sizeby == "null"){
    prob <- function(v) exp(
      lchoose(n, v[1]) + lchoose(n, v[2]) + dbetabinom(sum(v), 2*n, a, b, log = T) -
        lchoose(2*n, sum(v))
    )
  } else { # sizeby == "alt"
    prob <- function(v) exp(
      dbetabinom(v[1], n, a1, b1, log = T) + dbetabinom(v[2], n, a2, b2, log = T)
    )
  }

  df$p <- apply(df[,c("y1","y2")], 1, prob)
    
  qplot(y1, y2, data = df, geom = "point", color = action, shape = I(15), size = p) +
    scale_color_manual("Decision",
      values = c(
        "Accept" = rgb(61, 153, 86, maxColorValue = 256),
        "Reject" = rgb(255, 10, 27, maxColorValue = 256)
      )
    ) +
    scale_size("Probability")  
    
}



































#' Plot a Poisson decision rule
#'
#' Plot a Poisson decision rule
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
#' @param c1 loss associated with type I error
#' @param c2 loss associated with type I error
#' @param c relative loss constant (loss due to type II error divided by loss due to type I error)
#' @param rule the decision rule for the Bayes factor B.  by default, test rejects when B < pi1/pi0 * c2/c1
#' @param sizeby size the outcomes by the probability under the null or alternative
#' @param geom ggplot2 geom used
#' @return a ggplot object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @export plotPoisRule
#' @examples
#' 
#' \dontrun{
#'  
#' plotGamma(c(1,1), c(1,1))
#' plotPoisRule(1, 1, 1, 1, 1)
#' plotPoisRule(2, 1, 1, 1, 1)
#' plotPoisRule(3, 1, 1, 1, 1)
#' plotPoisRule(4, 1, 1, 1, 1)
#' plotPoisRule(5, 1, 1, 1, 1)
#' plotPoisRule(10, 1, 1, 1, 1)
#' plotPoisRule(5, 1, 1, 1, 1, sizeby = "null")
#' plotPoisRule(5, 1, 1, 1, 1, sizeby = "alt") # uniform
#' 
#' plotPoisRule(30, 3, 7, 7, 3)
#' plotPoisRule(100, 3, 7, 7, 3)
#'
#' plotPoisRule(30, 3, 7, 7, 3, sizeby = "null") + theme_bw()
#' plotPoisRule(30, 3, 7, 7, 3, sizeby = "alt") + theme_bw()
#'
#' }
#' 
plotPoisRule <- function(t, a1, b1, a2, b2, 
  a = a1, b = b1, pi0 = .5, pi1 = 1 - pi0, c1 = 1, c2 = 1, c = c1/c2, rule = pi1/pi0 * c2/c1,
  sizeby = c("none", "null", "alt"), geom = c("auto", "point", "tile"))
{
  
  # fake R CMD check
  y1 <- NULL; rm(y1);
  y2 <- NULL; rm(y2);  
  p <- NULL; rm(p);    
  action <- NULL; rm(action);  
  sizeby <- match.arg(sizeby)
  geom <- match.arg(geom)  
  
  if(sizeby != "none") geom <- "point"
  if(sizeby != "none") stop("this method is not yet implemented.", call. = FALSE)
  
  y1min <- qnbinom(.00001, a1, b1/(t+b1))
  y2min <- qnbinom(.00001, a2, b2/(t+b2))
  
  y1max <- qnbinom(.99999, a1, b1/(t+b1))
  y2max <- qnbinom(.99999, a2, b2/(t+b2))
  
  df <- expand.grid(y1 = y1min:y1max, y2 = y2min:y2max)
  
  if(geom == "auto" && nrow(df) <= 50^2){
    geom <- "point"
  } else if(geom == "auto" && nrow(df) > 50^2){
    geom <- "tile"    
  }
  
  df$action <- apply(df, 1, function(x) 
    bayesPoisTest(x, t, a1, b1, a2, b2, a, b, pi0, pi1, c, rule)$reject
  )
  
  df$action[df$action == TRUE] <- "Reject"
  df$action[df$action == "FALSE"] <- "Accept"  
  
  if(geom == "tile") return(
    qplot(y1, y2, data = df, geom = "tile", fill = action) +
      scale_fill_manual("Decision",
        values = c(
          "Accept" = rgb(61, 153, 86, maxColorValue = 256),
          "Reject" = rgb(255, 10, 27, maxColorValue = 256)
        )
      )
  )
  
  if(sizeby == "none") return(    
    qplot(y1, y2, data = df, geom = "point", color = action, shape = I(15)) +
      scale_color_manual("Decision",
        values = c(
          "Accept" = rgb(61, 153, 86, maxColorValue = 256),
          "Reject" = rgb(255, 10, 27, maxColorValue = 256)
        )
      )    
  )
  
  n <- NULL; rm(n);    
  if(sizeby == "null"){
    prob <- function(v) exp(
      lchoose(n, v[1]) + lchoose(n, v[2]) + dbetabinom(sum(v), 2*n, a, b, log = T) -
        lchoose(2*n, sum(v))
    )
  } else { # sizeby == "alt"
    prob <- function(v) exp(
      dbetabinom(v[1], n, a1, b1, log = T) + dbetabinom(v[2], n, a2, b2, log = T)
    )
  }
  
  df$p <- apply(df[,c("y1","y2")], 1, prob)
  
  qplot(y1, y2, data = df, geom = "point", color = action, shape = I(15), size = p) +
    scale_color_manual("Decision",
      values = c(
        "Accept" = rgb(61, 153, 86, maxColorValue = 256),
        "Reject" = rgb(255, 10, 27, maxColorValue = 256)
      )
    ) +
    scale_size("Probability")  
  
}


