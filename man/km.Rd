\name{km}
\title{Kaplan-Meier Survival Curves}
\usage{
km(times, censor, group=1, freq=1, cdf=F)
}
\alias{km}
\alias{plot.hazard}
\alias{plot.dist}
\arguments{
\item{times}{Vector of times to events or a list of vectors of such
times for different individuals.}
\item{censor}{Vector of censoring indicators corresponding to the
vector of times or to the last time in each vector of a list.}
\item{group}{Vector indicating to which group each individual belongs.}
\item{freq}{Vector of frequencies for grouped data.}
\item{cdf}{If TRUE, calculate the cdf instead of the survival curve.}}
}
\description{
\code{km} calculates the Kaplan-Meier estimates for survival.

To plot the survival curve, use plot(); for the empirical hazard curve,
use plot.hazard(); for diagnostic curves to choose a distribution to
which the data might belong, use plot.dist().
}
\value{
A matrix with class, km, containing the Kaplan-Meier estimates
is returned.
}
\author{J.K. Lindsey}
\examples{
surv <- rgamma(40,2,5)
cens <- rbinom(40,1,0.9)
treat <- gl(2,20)
plot(km(surv, cens, group=treat), main="",xlab="Months",
	ylab="Probability of deterioration")
plot.dist(km(surv, cens, group=treat))
plot.hazard(km(surv, cens, group=treat),ylab="Risk of deterioration")
}
