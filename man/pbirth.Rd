\name{pbirth}
\title{Fit Overdispersed Count Data as a Birth Process}
\usage{
pbirth(frequencies, p, intensity="negative binomial",
	type="spectral decomposition", print.level=0, typsiz=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p\%*\%p),
	steptol=0.00001, iterlim=100, fscale=1)
}
\alias{pbirth}
\arguments{
\item{frequencies}{Vector of frequencies or a matrix with each row a
different series of frequencies.}
\item{p}{Vector of initial estimates.}
\item{intensity}{The intensity function of the process: binomial,
binomial exdponential, binomial logistic, binomial total,
Poisson, Poisson exponential, negative binomial, or gen(eralized)
negative binomial.}
\item{type}{Algorithm used for matrix exponentiation: spectral
decomposition or series approximation.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{pbirth} fits binomial, binomial exponential, binomial logistic,
binomial total, Poisson, Poisson exponential, negative binomial,
gen(eralized) negative binomial, and generalized negative binomial
processes as a birth process.
}
\references{
Faddy, M.J. and Fenlon, J.S. (1999) Stochastic modelling of the
invasion process of nematodes in fly larvae. Applied Statistics 48:
31-37.
}
\author{J.K. Lindsey}
\examples{
y <- rnbinom(100,2,0.6)
fr <- tabulate(y)
pbirth(fr, p=log(-log(0.7)), intensity="Poisson", type="series")
pbirth(fr, p=c(log(-log(0.7)),log(5)),
	intensity="negative binomial", type="series")
pbirth(fr, p=c(log(-log(0.7)),log(5),-1),
	intensity="gen negative binomial", type="series")
}
