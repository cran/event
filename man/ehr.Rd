\name{ehr}
\title{Fit an Intensity Function to Event Histories}
\usage{
ehr(point, lambda=NULL, linear=NULL, plambda=NULL, delta=1,
	envir=sys.frame(sys.parent()), print.level=0,
	typsiz=rep(1,length(p)), ndigit=10, gradtol=0.00001,
	stepmax=max(10*sqrt(p\%*\%p),10), steptl=0.0004,
	iterlim=100, fscale=1)
}
\alias{ehr}
\arguments{
\item{point}{A point process vector produced by \code{\link{pp}}.}
\item{lambda}{User-specified function of \code{p}, and possibly
\code{linear}, giving the regression equation for the intensity or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation or a general function
with named unknown parameters. The function may contain a
linear part that must simply be given the name, \code{linear}, in the
function. If no function is supplied, the intensity is taken to be
constant (a homogeneous Poisson process).}
\item{linear}{A formula beginning with ~ specifying the linear part of
the regression function.}
\item{plambda}{Vector of initial parameter estimates. If \code{lambda}
is a formula with unknown parameters, their estimates must be supplied
either in their order of appearance in the expression or in a named list.}
\item{delta}{If any time intervals are different from unity, a vector
of time intervals.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{ehr} fits an intensity function to event histories, where point is
produced by \code{point <- pp(y)} and \code{lambda} is the user-defined
log intensity function.

Nonlinear regression models for \code{lambda} can be supplied as
formulae where parameters are unknowns. Factor variables cannot be
used and parameters must be scalars. (See \code{\link{finterp}}.)
}
\references{
Lindsey, J.K. (1995) Fitting parametric counting processes by
using log linear models. Journal of the Royal Statistical
Society C44, 201-212.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{finterp}}, \code{\link{ident}},
\code{\link{pp}}, \code{\link{tccov}}, \code{\link{tpast}},
\code{\link{ttime}}, \code{\link{tvcov}}.
}
\examples{
y <- c(5,3,2,4)
# event indicator
py <- pp(y)
# time since previous event
ptime <- tpast(y)
# individual ID
i <- c(1,1,2,2)
id <- ident(y, i)
# times and corresponding covariate values
tx <- c(2,3,1,2,2,2,2)
x <- c(1,2,2,1,2,2,1)
zcov <- tvcov(y, x, tx)
# Poisson process
ehr(py, plambda=1)
# Weibull process
lambda1 <- function(p) p[1]+p[2]*log(ptime)
ehr(py, lambda=lambda1, plambda=c(1,1))
# or
ehr(py, lambda=~log(ptime), plambda=c(1,1))
# or
ehr(py, lambda=~b0+b1*log(ptime), plambda=list(b0=1,b1=1))
# Poisson process with time-varying covariate
lambda2 <- function(p) p[1]+p[2]*zcov
ehr(py, lambda=lambda2, plambda=c(1,1))
# or
ehr(py, lambda=~zcov, plambda=c(1,1))
# or
ehr(py, lambda=~c0+c1*zcov, plambda=list(c0=1,c1=1))
# Weibull process with time-varying covariate
lambda3 <- function(p) p[1]+p[2]*log(ptime)+p[3]*zcov
ehr(py, lambda=lambda3, plambda=c(1,1,1))
# or
ehr(py, lambda=~log(ptime)+zcov, plambda=c(1,1,1))
# or
ehr(py, lambda=~c0+b1*log(ptime)+c1*zcov, plambda=list(c0=1,c1=1,b1=1))
# gamma process with time-varying covariate
lambda4 <- function(p) hgamma(ptime, p[1], exp(p[2]+p[3]*zcov))
ehr(py, lambda=lambda4, plambda=c(1,1,1))
# or
ehr(py, lambda=~hgamma(ptime, b1, exp(c0+c1*zcov)),
	plambda=list(c0=1,c1=1,b1=1))
}
