\name{pp}
\title{Create a Point Process Vector from Times between Events}
\usage{
pp(y, censor=1)
}
\alias{pp}
\arguments{
\item{y}{Vector of times.}
\item{cens}{Vector of censoring indicators.}
}
\value{
\code{pp} creates a vector of length \code{sum(y)} of zeroes with a one
at the end of each uncensored time interval for use with \code{\link{ehr}}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{ehr}}, \code{\link{ident}}, \code{\link{tccov}},
\code{\link{tpast}}, \code{\link{ttime}}, \code{\link{tvcov}}.
}
\examples{
y <- c(5,3,2,4)
py <- pp(y)
py
}
