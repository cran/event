\name{tvcov}
\title{Create a Vector of Time-varying Covariates for
a Point Process}
\usage{
tvcov(y, x, tx)
}
\alias{tvcov}
\arguments{
\item{y}{Vector of times.}
\item{x}{Vector covariate.}
\item{tx}{Vector of times at which x changes.}
}
\value{
\code{tvcov} creates a vector of length \code{sum(y)} of time-varying
covariates for use with \code{\link{ehr}}. It may also be used to create
weight vectors.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{ehr}}, \code{\link{ident}}, \code{\link{pp}},
\code{\link{tccov}}, \code{\link{tpast}}, \code{\link{ttime}}.
}
\examples{
y <- c(5,3,2,4)
x <- c(1,2,2,1,2,2,1)
tx <- c(2,3,1,2,2,2,2)
zcov <- tvcov(y, x, tx)
zcov
}
