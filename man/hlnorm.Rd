\name{hlnorm}
\title{Hazard Function for a Log Normal Process}
\usage{
hlnorm(y, m, s)
}
\alias{hlnorm}
\arguments{
\item{y}{Vector of times.}
\item{m}{Mean parameter.}
\item{s}{Variance parameter.}
}
\value{
\code{hlnorm} returns the log hazard function for a log normal
process with the given parameter values.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{hboxcox}}, \code{\link{hburr}}, \code{\link{hexp}},
\code{\link{hgextval}}, \code{\link{hcauchy}}, \code{\link{hgamma}},
\code{\link{hggamma}}, \code{\link{hhjorth}}, \code{\link{hinvgauss}},
\code{\link{hlaplace}}, \code{\link{hlogis}}, \code{\link{hglogis}},
\code{\link{hnorm}}, \code{\link{hstudent}}, \code{\link{hweibull}},
\code{\link{hgweibull}}.
}
\examples{
hlnorm(1:10, 3, 2)
}
