\name{ident}
\title{Create an Individual Identification Vector for a Point Process}
\usage{
ident(y, id)
}
\alias{ident}
\arguments{
\item{y}{Vector of times.}
\item{id}{Vector of corresponding individual identifiers for who had
which sequence of times.}
}
\value{
\code{ident} creates a vector of length \code{sum(y)} by repeating the
values of individual identifiers for the times for use with \code{\link{ehr}}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{bp}}, \code{\link{ehr}}, \code{\link{pp}}, \code{\link{tccov}},
\code{\link{tpast}}, \code{\link{ttime}}, \code{\link{tvcov}}.
}
\examples{
y <- c(5,3,2,4)
i <- c(1,1,2,2)
id <- ident(y, i)
id
}
