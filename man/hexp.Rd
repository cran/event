\name{hexp}
\title{Hazard Function for a Poisson Process}
\usage{
hexp(y, m)
}
\alias{hexp}
\arguments{
\item{y}{Vector of times.}
\item{m}{Mean parameter.}
}
\value{
\code{hexp} returns the log hazard function for a Poisson process with
the given parameter value.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{hboxcox}}, \code{\link{hburr}}, \code{\link{hexp}},
\code{\link{hgextval}}, \code{\link{hgamma}}, \code{\link{hggamma}},
\code{\link{hhjorth}}, \code{\link{hinvgauss}}, \code{\link{hlaplace}},
\code{\link{hlnorm}}, \code{\link{hlogis}}, \code{\link{hglogis}},
\code{\link{hnorm}}, \code{\link{hstudent}}, \code{\link{hweibull}},
\code{\link{hgweibull}}.
}
\examples{
hexp(1:10, 3)
}
