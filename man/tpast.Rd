\name{tpast}
\title{Create a Vector of Times Past since Previous Events for
a Point Process}
\usage{
tpast(y)
}
\alias{tpast}
\arguments{
\item{y}{Vector of times.}
}
\value{
\code{tpast} creates a vector of length \code{sum(y)} of times since the
immediately preceding event occurred for use with \code{\link{ehr}}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{ehr}}, \code{\link{ident}}, \code{\link{pp}},
\code{\link{tccov}}, \code{\link{ttime}}, \code{\link{tvcov}}.
}
\examples{
y <- c(5,3,2,4)
ptime <- tpast(y)
ptime
}
