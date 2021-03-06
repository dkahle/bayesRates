% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.r
\name{bayesPoisTest}
\alias{bayesPoisTest}
\title{Bayesian two-sample test of Poisson rates}
\usage{
bayesPoisTest(x, t, a1, b1, a2, b2, a = a1, b = b1, pi0 = 0.5,
  pi1 = 1 - pi0, c = 1, rule = c * pi1/pi0)
}
\arguments{
\item{x}{two-vector containing occurrences}

\item{t}{duration of trials (the same for both samples)}

\item{a1}{alpha, the hyperparameter of the gamma distribution of the first poisson rate}

\item{b1}{beta, the hyperparameter of the gamma distribution of the first poisson rate}

\item{a2}{alpha, the hyperparameter of the gamma distribution of the second poisson rate}

\item{b2}{beta, the hyperparameter of the gamma distribution of the second poisson rate}

\item{a}{alpha, the hyperparameter of the gamma distribution under the null}

\item{b}{beta, the hyperparameter of the gamma distribution under the null}

\item{pi0}{the prior probability of the null hypothesis}

\item{pi1}{the prior probability of the alternative hypothesis}

\item{c}{relative loss constant (loss due to type II error divided by loss due to type I error)}

\item{rule}{the decision rule for the Bayes factor B.  by default, test rejects when B < pi1/pi0 * c2/c1}
}
\description{
Bayesian two-sample test of Poisson rates
}
\examples{
\dontrun{

a1 <- 5; b1 <- .2 # E = 1
a2 <- 30; b2 <- .167 # E = 5
plotGamma(c(a1, a2), c(b1, b2))
bayesPoisTest(c(18, 105), 20, a1, b1, a2, b2)
bayesPoisTest(c(18, 165), 20, a1, b1, a2, b2)

# compare to :
prop.test(c(10, 16), c(20, 20))


 }
}
\seealso{
\code{\link{samplePower}}, \code{\link{samplePowerEst}}, \code{\link{findSize}}
}
