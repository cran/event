\name{hweibull}
\title{Hazard Function for a Weibull Process}
\usage{
hweibull(y, s, m)
}
\alias{hweibull}
\arguments{
\item{y}{Vector of times.}
\item{s}{Shape parameter.}
\item{m}{Scale parameter.}
}
\value{
\code{hweibull} returns the log hazard function for a Weibull
process with the given parameter values.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{hboxcox}}, \code{\link{hburr}}, \code{\link{hexp}},
\code{\link{hgextval}}, \code{\link{hcauchy}}, \code{\link{hgamma}},
\code{\link{hggamma}}, \code{\link{hhjorth}}, \code{\link{hinvgauss}},
\code{\link{hlaplace}}, \code{\link{hlogis}}, \code{\link{hglogis}},
\code{\link{hnorm}}, \code{\link{hnorm}}, \code{\link{hstudent}},
\code{\link{hgweibull}}.
}
\examples{
hweibull(1:10, 1.5, 2)
}
