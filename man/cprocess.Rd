\name{cprocess}
\title{Counting Process Plot}
\usage{
cprocess(times=NULL, events=NULL, number=F, lty=NULL, xlim=NULL,
	ylim=NULL, xlab="Time", ylab="Counts", ...)
}
\alias{cprocess}
\arguments{
\item{times}{Vector of times between events or a list of such
vectors.}
\item{events}{Vector of counts of events or a list of such vectors.}
\item{number}{If TRUE, the counting processes are numbered
consecutively on the graph.}
\item{others}{Plotting control options.}
}
\description{
\code{cprocess} plots the cumulative number of events (the counting
process) over time for each vector in the list. One or both of
\code{times} and \code{events} can be supplied.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{read.list}}.
}
\examples{
times <- rgamma(20,2,4)
cprocess(times)
}
