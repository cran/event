\name{coxre}
\title{Cox Proportional Hazards Model with Random Effect}
\usage{
coxre(repsonse, censor, nest, cov, stratified=F, cumul=F, estimate=1,
	iter=10, print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, iterlim=100, fscale=1, typsiz=abs(estimate),
	stepmax=estimate)
}
\alias{coxre}
\arguments{
\item{response}{Vector or matrix of times to events, with one column
per type of response (or subunit).}
\item{censor}{Corresponding vector or matrix of censoring indicators.}
\item{nest}{Vector indicating to which unit each observation belongs.}
\item{cov}{One covariate}
\item{stratified}{If TRUE, a model stratified on type of response
(the columns of response) is fitted instead of proportional intensities.}
\item{cumul}{Set to TRUE if response times are from a common origin
instead of times to (or between) events.}
\item{estimate}{Initial estimate of the frailty parameter.}
\item{iter}{Maximum number of iterations allowed.}
\item{others}{Plotting control options.}
}
\description{
\code{coxre} fits a Cox proportional hazards model to event history
data using a gamma distribution random effect. The parameter, gamma,
is the variance of this mixing distribution.

If a matrix of response times is supplied, the model can be stratified
by columns, i.e. a different intensity function is fitted for each
column. To fit identical intensity functions to all response types,
give the times as a vector.
}
\author{D.G. Clayton and J.K. Lindsey}
\references{
Clayton, D. (1987) The analysis of event history data: a review of
progress and outstanding problems.  Statistics in Medicine 7: 819-841
}
\seealso{
\code{\link{kalsurv}}.
}
\examples{
# 11 individuals, each with 5 responses
y <- matrix(rweibull(55,2,5),ncol=5)
# Different intensity functions
coxre(response=y, censor=matrix(rep(1,55),ncol=5), nest=1:11,
	est=0.7, stratified=T)
# Proportional intensity functions
coxre(response=y, censor=matrix(rep(1,55),ncol=5),
	nest=1:11, cov=rpois(11,2), est=0.7, stratified=F)
# Identical intensity functions
coxre(response=as.vector(t(y)), censor=rep(1,55),
	nest=rep(1:11,rep(5,11)), est=0.7)
}
