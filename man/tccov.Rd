\name{tccov}
\title{Create a Vector of Time-constant Covariates for
a Point Process}
\usage{
tccov(y, x, id)
}
\alias{tccov}
\arguments{
\item{y}{Vector of times.}
\item{x}{Vector covariate.}
\item{id}{Vector of corresponding individual identifiers for who had
which sequence of times.}
}
\value{
\code{tccov} creates a vector of length \code{sum(y)} of time-constant
covariates for use with \code{\link{ehr}}. \code{id} must be numbered
consecutively. \code{x} must have one value for each distinct \code{id},
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{ehr}}, \code{\link{ident}}, \code{\link{pp}},
\code{\link{tpast}}, \code{\link{ttime}}, \code{\link{tvcov}}.
}
\examples{
y <- c(5,3,2,4)
id <- c(1,1,2,2)
x <- c(5.2,3.1)
xcov <- tccov(y, x, id)
xcov
}
