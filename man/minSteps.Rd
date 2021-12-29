\name{minSteps}
\alias{minSteps}
\title{minSteps}
\description{Calculates the minimum amount of times must execute bissection method to get some "good"" value.}
\usage{minSteps(a,b,E)}
\arguments{
  \item{a}{\strong{INT}
        Lower bound}
  \item{b}{\strong{INT}
        Upper bound}
  \item{E}{\strong{INT}
        Tolerance}
}
\value{\strong{INT} ceiling(value)}
\author{@ppcamp}
\note{Classnotes by João}
\examples{minSteps(a=1,b=2,E=10^-3)}
\keyword{ root }% __ONLY ONE__ keyword per line
