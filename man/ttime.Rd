\name{ttime}
\title{Create a Vector of Total Time Elapsed for each
Individual for a Point Process}
\usage{
ttime(y, id)
}
\alias{ttime}
\arguments{
\item{y}{Vector of times.}
\item{id}{Vector of corresponding individual identifiers for who had
which sequence of times.}
}
\value{
\code{ttime} creates a vector of length \code{sum(y)} of times since
each individual began for use with \code{\link{ehr}}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{ehr}}, \code{\link{ident}}, \code{\link{pp}},
\code{\link{tccov}}, \code{\link{tpast}}, \code{\link{tvcov}}.
}
\examples{
y <- c(5,3,2,4)
id <- c(1,1,2,2)
itime <- ttime(y, id)
itime
}
