\name{hgamma}
\title{Hazard Function for a Gamma Process}
\usage{
hgamma(y, s, m)
}
\alias{hgamma}
\arguments{
\item{y}{Vector of times.}
\item{s}{Shape parameter.}
\item{m}{Scale parameter.}
}
\value{
\code{hgamma} returns the log hazard function for a gamma process with
the given parameter values.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{hboxcox}}, \code{\link{hburr}}, \code{\link{hexp}},
\code{\link{hgextval}}, \code{\link{hcauchy}}, \code{\link{hggamma}},
\code{\link{hhjorth}}, \code{\link{hinvgauss}}, \code{\link{hlaplace}},
\code{\link{hlnorm}}, \code{\link{hlogis}}, \code{\link{hglogis}},
\code{\link{hnorm}}, \code{\link{hstudent}}, \code{\link{hweibull}},
\code{\link{hgweibull}}.
}
\examples{
hgamma(1:10, 3, 2)
}
