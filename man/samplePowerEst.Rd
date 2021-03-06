% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplePowerEst.r
\name{samplePowerEst}
\alias{samplePowerEst}
\title{Approximate the power of a given sample size}
\usage{
samplePowerEst(n, a1, b1, a2, b2, a = a1, b = b1, pi0 = 0.5,
  pi1 = 1 - pi0, method = c("fastest", "fast", "moderate"),
  family = c("binomial", "poisson"))
}
\arguments{
\item{n}{the sample size - possibly a vector}

\item{a1}{alpha, the hyperparameter of the gamma distribution of the first poisson rate}

\item{b1}{beta, the hyperparameter of the gamma distribution of the first poisson rate}

\item{a2}{alpha, the hyperparameter of the gamma distribution of the second poisson rate}

\item{b2}{beta, the hyperparameter of the gamma distribution of the second poisson rate}

\item{a}{alpha, the hyperparameter of the gamma distribution under the null}

\item{b}{beta, the hyperparameter of the gamma distribution under the null}

\item{pi0}{the prior probability of the null hypothesis}

\item{pi1}{the prior probability of the alternative hypothesis}

\item{method}{"fastest", "fast", "moderate", decreasing in speed and increasing in accuracy}

\item{family}{"binomial" or "poisson"}
}
\description{
Approximate a sample size
}
\examples{
\dontrun{

# binomial
a1 <- 2; b1 <- 6
a2 <- 3; b2 <- 3
samplePower(n = 10, a1, b1, a2, b2)
samplePowerEst(n = 10, a1, b1, a2, b2)

# poisson
a1 <- 4; b1 <- 2
a2 <- 6; b2 <- 2
samplePower(n= 20, a1, b1, a2, b2, family = "poisson")
samplePowerEst(n = 20, a1, b1, a2, b2, family = "poisson")
samplePowerEst(n = 20, a1, b1, a2, b2, family = "poisson", method = "fast")
samplePowerEst(n = 20, a1, b1, a2, b2, family = "poisson", method = "moderate")

# this takes a moment
samplePowerEst(n = 100*1:10, a1, b1, a2, b2, family = "poisson")
samplePower(n = 200, a1, b1, a2, b2, family = "poisson")









	
## predict power for a given sample size (poisson)
############################################################
n <- 100
mult <- 1
power <- numeric(n)
time <- numeric(n)
t <- c(mult*(1:n),150,200,250,300,350,400,450,500,750,1000)
t <- c(1, 50, 100 ,150,200,250,300,350,400,450,500,750,1000)
for(k in seq_along(t)){
  start <- Sys.time()
  message(paste0("t = ", t[k], "... "), appendLF = FALSE)   
  power[k] <- samplePower(t[k], a1, b1, a2, b2, family = "poisson")
  time[k] <- end <- difftime(Sys.time(), start, units = "secs")
  if(k == length(t)) message("")
}
df <- data.frame(n = t, time = time, power = power)


library(ggplot2)
# investigate power as a function of sample size n
qplot(n, power, data = df) + ylim(0,1)

qplot(n, power, data = df) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  ylim(0,1)
  
mod <- lm(power ~ log(n), data = df)  
s <- 1:1000
pdf <- data.frame(
  n = s,
  power = predict(mod, data.frame(n = s))
)  
qplot(n, power, data = pdf, geom = "path", ylim = c(0,1)) +
    labs(x = "n (Sample Size)", y = "Estimated Power") +
  geom_point(aes(x = n, y = power), data = df, color = "red")



## predict time necessary for the exact power to be computed as a 
## function of the sample size (poisson)
############################################################

qplot(n, time, data = df)

# not log...
qplot(n, time, data = df) + scale_y_log10()

# quadratic
qplot(n, time, data = df) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2))
  
mod <- lm(time ~ n + I(n^2), data = df)
s <- 1:1000
tdf <- data.frame(
  n = s,
  time = predict(mod, data.frame(n = s))
)  
qplot(n, time, data = tdf, geom = "path") +
  labs(x = "n (Sample Size)", y = "Time (Seconds)")


start <- Sys.time()
samplePower(500, a1, b1, a2, b2, family = "poisson")
difftime(Sys.time(), start, units = "secs")
subset(tdf, n == 500)
# predicted times are very close



}
}
\seealso{
\code{\link{samplePower}}, \code{\link{findSize}}
}
