\name{hcauchy}
\title{Hazard Function for a Cauchy Process}
\usage{
hcauchy(y, m, s)
}
\alias{hcauchy}
\arguments{
\item{y}{Vector of times.}
\item{m}{Location parameter.}
\item{s}{Dispersion parameter.}
}
\value{
\code{hcauchy} returns the log hazard function for a Cauchy process with
the given parameter values.
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
hcauchy(1:10, 3, 2)
}
