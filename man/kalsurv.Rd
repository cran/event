\name{kalsurv}
\title{Repeated Events Models for Changes over Time in the
Presence of Time-varying Covariates and with Frailty or Various Types
of Serial Dependence}
\usage{
kalsurv(response, intensity="exponential", distribution="Pareto",
	depend="independence", update="Markov", mu=NULL, shape=NULL,
	renewal=T, density=F, censor=NULL, delta=NULL, ccov=NULL,
	tvcov=NULL, preg=1, ptvc=NULL, pbirth=NULL,
	pintercept=NULL, pshape=NULL, pinitial=1, pdepend=NULL,
	envir=sys.frame(sys.parent()), print.level=0,
	ndigit=10, gradtol=0.00001, steptol=0.00001, iterlim=100,
	fscale=1, typsiz=abs(p), stepmax=10*sqrt(p\%*\%p))
}
\alias{kalsurv}
\alias{coefficients.kalsurv}
\alias{deviance.kalsurv}
\alias{fitted.kalsurv}
\alias{residuals.kalsurv}
\alias{print.kalsurv}
\arguments{
\item{response}{A list of vectors with times between events for
each individual, one matrix or dataframe of such times if all
individuals have the same number of events, or an object of class,
response (created by \code{\link{restovec}}).}
\item{intensity}{The form of intensity function to be put in the
distribution given by dist. Choices are exponential, Weibull, gamma,
log normal, log logistic, log Cauchy, log Student, and gen(eralized)
logistic.}
\item{distribution}{The outer distribution. Choices are Pareto, gamma,
and Weibull.}
\item{depend}{Type of dependence. Choices are independence, frailty,
and serial.}
\item{update}{Type of update for serial dependence. Choices are
Markov, elapsed Markov, serial, event, cumulated, count, and
kalman. With frailty dependence, weighting by length of observation
time may be specified by setting update to \code{time}.}
\item{mu}{A regression function for the location parameter or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation or a general function
with named unknown parameters. If there are only time-constant
covariates, give the initial estimates in preg; if any covariates are
time-varying, give all initial estimates in ptvc.}
\item{shape}{A regression function for the shape parameter or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation or a general function
with named unknown parameters. It must yield one value per observation.}
\item{renewal}{IF TRUE, a renewal process is modelled, with time
reinitialized after each event. Otherwise, time is cumulated from the
origin of observations.}
\item{density}{If TRUE, the density of the function specified in
\code{intensity} is used instead of the intensity.}
\item{censor}{A vector of the same length as the number of individuals
containing a binary indicator, with a one indicating that the last
time period in the series terminated with an event and zero that it
was censored. For independence and frailty models, where response is
matrix, censor may also be a matrix of the same size. Ignored if
response has class, response or repeated.}
\item{delta}{Scalar or vector giving the unit of measurement for each
response value, set to unity by default. For example, if a response is
measured to two decimals, delta=0.01. If the response has been
pretransformed, this must be multiplied by the Jacobian. This
transformation cannot contain unknown parameters. For example, with a
log transformation, \code{delta=1/y}. (The delta values for the
censored response are ignored.) Ignored if response has class,
response or repeated.}
\item{ccov}{A vector or matrix containing time-constant baseline
covariates with one entry per individual, a model formula using
vectors of the same size, or an object of class, tccov (created by
\code{\link{tcctomat}}). Ignored if response has class, repeated.}
\item{tvcov}{A list of matrices with time-varying covariate values,
observed at the event times in \code{response}, for each individual
(one column per variable), one matrix or dataframe of such covariate
values, or an object of class, tvcov (created by
\code{\link{tvctomat}}). Ignored if response has class, repeated.}
\item{preg}{Initial parameter estimates for the regression model:
intercept plus one for each covariate in \code{ccov}. If a location
function (mu) is supplied that contains time-varying covariates, all
initial estimates must be given in ptvc.}
\item{ptvc}{Initial parameter estimates for the coefficients of the
time-varying covariates, as many as in \code{tvcov}. If a location
function (mu) is supplied that contains time-varying covariates, all
initial estimates must be given here.}
\item{pbirth}{If supplied, this is the initial estimate for the
coefficient of the birth model.}
\item{pintercept}{The initial estimate of the intercept for the
generalized logistic intensity.}
\item{pshape}{An initial estimate for the shape parameter of the
intensity (except exponential intensity). If \code{shape} is
a function or formula, the corresponding initial estimates. If
\code{shape} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pinitial}{An initial estimate for the initial parameter. (In
frailty dependence, this is the frailty parameter.)}
\item{pdepend}{An initial estimate for the serial dependence
parameter.}
\item{envir}{Environment in which model formulae are to be
interpreted.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{kalsurv} is designed to handle event history models with
time-varying covariates. The distributions have two extra parameters
as compared to the functions specified by \code{intensity} and are
generally longer tailed than those distributions. Dependence of
inter-event times can be through frailty (a type of random effect) or
several kinds of serial dependence by updating, as in Kalman
filtering.

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)

Marginal and individual profiles can be plotted using
\code{\link{plot.profile}} and \code{\link{plot.iprofile}} and
residuals with \code{\link{plot.residuals}}.
}
\value{
A list of classes \code{kalsurv} and \code{recursive} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{contkal}}, \code{\link{coxre}}, \code{\link{finterp}},
\code{\link{gettvc}}, \code{\link{gnlmm}}, \code{\link{gnlr}},
\code{\link{kalcount}}, \code{\link{kalseries}}, \code{\link{nbkal}},
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{rmna}},
\code{\link{tcctomat}}, \code{\link{tvctomat}}.
}
\examples{
treat <- c(0,0,1,1)
tr <- tcctomat(treat)
cens <- matrix(rbinom(20,1,0.9),ncol=5)
times <- # matrix(rweibull(20,2,1+3*rep(treat,5)),ncol=5)
	matrix(c(1.36,0.18,0.84,0.65,1.44,1.79,1.04,0.43,1.35,1.63,2.15,1.15,
		1.21,5.46,1.58,3.44,4.40,2.75,4.78,2.44),ncol=5,byrow=T)
times <- restovec(times, censor=cens)
reps <- rmna(times, ccov=tr)
# exponential intensity model with independence
kalsurv(times, pinitial=0.1, preg=1,
	dep="independence", intensity="exponential")
# Weibull intensity model with independence
kalsurv(times, pinitial=0.1, preg=1, pshape=1,
	dep="independence", intensity="Weibull")
# same model with serial update
kalsurv(times, pinitial=0.1, pdep=0.1, preg=1, pshape=1,
	dep="serial", intensity="Weibull")
# treatment effect with log link
kalsurv(times, pinitial=0.1, preg=c(1,0), pshape=1,
	intensity="Weibull", ccov=treat)
# or equivalently
kalsurv(times, mu=~exp(a+b*treat), pinitial=0.1, preg=c(1,0), pshape=1,
	intensity="Weibull", envir=reps)
# with identity link instead
kalsurv(times, mu=~treat, pinitial=0.1, preg=c(1,0), pshape=1,
	intensity="Weibull")
# or equivalently
kalsurv(times, mu=~a+b*treat, pinitial=0.1, preg=c(1,0), pshape=1,
	intensity="Weibull", envir=reps)
# add the birth model
kalsurv(times, pinitial=0.1, preg=c(1,0), pshape=1,
	intensity="Weibull", ccov=treat, pbirth=0)
# switch to gamma intensity model without birth effect
kalsurv(times, pinitial=0.1, preg=c(1,0), pshape=1, intensity="gamma",
	ccov=treat)
# try frailty dependence
kalsurv(times, pinitial=0.1, preg=c(1,0), pshape=1, dep="frailty",
	intensity="gamma", ccov=treat)
}
