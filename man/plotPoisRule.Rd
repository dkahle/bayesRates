% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotRule.r
\name{plotPoisRule}
\alias{plotPoisRule}
\title{Plot a Poisson decision rule}
\usage{
plotPoisRule(t, a1, b1, a2, b2, a = a1, b = b1, pi0 = 0.5, pi1 = 1
  - pi0, c1 = 1, c2 = 1, c = c1/c2, rule = pi1/pi0 * c2/c1,
  sizeby = c("none", "null", "alt"), geom = c("auto", "point", "tile"))
}
\arguments{
\item{t}{the sample size}

\item{a1}{alpha, the hyperparameter of the gamma distribution of the first poisson rate}

\item{b1}{beta, the hyperparameter of the gamma distribution of the first poisson rate}

\item{a2}{alpha, the hyperparameter of the gamma distribution of the second poisson rate}

\item{b2}{beta, the hyperparameter of the gamma distribution of the second poisson rate}

\item{a}{alpha, the hyperparameter of the gamma distribution under the null}

\item{b}{beta, the hyperparameter of the gamma distribution under the null}

\item{pi0}{the prior probability of the null hypothesis}

\item{pi1}{the prior probability of the alternative hypothesis}

\item{c1}{loss associated with type I error}

\item{c2}{loss associated with type I error}

\item{c}{relative loss constant (loss due to type II error divided by loss due to type I error)}

\item{rule}{the decision rule for the Bayes factor B.  by default, test rejects when B < pi1/pi0 * c2/c1}

\item{sizeby}{size the outcomes by the probability under the null or alternative}

\item{geom}{ggplot2 geom used}
}
\value{
a ggplot object
}
\description{
Plot a Poisson decision rule
}
\examples{

\dontrun{
 
plotGamma(c(1,1), c(1,1))
plotPoisRule(1, 1, 1, 1, 1)
plotPoisRule(2, 1, 1, 1, 1)
plotPoisRule(3, 1, 1, 1, 1)
plotPoisRule(4, 1, 1, 1, 1)
plotPoisRule(5, 1, 1, 1, 1)
plotPoisRule(10, 1, 1, 1, 1)
plotPoisRule(30, 1, 1, 1, 1)
plotPoisRule(100, 1, 1, 1, 1)


plotPoisRule(5, 1, 1, 1, 1, sizeby = "null")
plotPoisRule(5, 1, 1, 1, 1, sizeby = "alt") 

plotPoisRule(30, 3, 7, 7, 3, geom = "point")
plotPoisRule(100, 3, 7, 7, 3)

plotPoisRule(30, 3, 7, 7, 3, sizeby = "null") + theme_bw()
plotPoisRule(30, 3, 7, 7, 3, sizeby = "alt") + theme_bw()

}

}
\author{
David Kahle \email{david.kahle@gmail.com}
}
