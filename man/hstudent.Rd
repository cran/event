\name{hstudent}
\title{Hazard Function for a Student t Process}
\usage{
hstudent(y, m, s, f)
}
\alias{hstudent}
\arguments{
\item{y}{Vector of times.}
\item{m}{Location parameter.}
\item{m}{Scale  parameter.}
\item{s}{Dispersion parameter.}
}
\value{
\code{hstudent} returns the log hazard function for a Student t
process with the given parameter values.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{hboxcox}}, \code{\link{hburr}}, \code{\link{hexp}},
\code{\link{hgextval}}, \code{\link{hcauchy}}, \code{\link{hgamma}},
\code{\link{hggamma}}, \code{\link{hhjorth}}, \code{\link{hinvgauss}},
\code{\link{hlaplace}}, \code{\link{hlogis}}, \code{\link{hglogis}},
\code{\link{hnorm}}, \code{\link{hnorm}}, \code{\link{hweibull}},
\code{\link{hgweibull}}.
}
\examples{
hstudent(1:10, 3, 2)
}
