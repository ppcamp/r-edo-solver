\name{fixedPointIteration}
\alias{fixedPointIteration}
\title{fixedPointIteration}
\description{Calculates root using fixed point iteration method.}
\arguments{
  \item{func}{\strong{STRING}
        Derivative formula write in STRING.
        NOTE: observes that exponential values must use
        '**' instead operator '^'}
  \item{p0}{\strong{NUMBER}
        Initial approximation}
  \item{tol}{\emph{OPTIONAL} \strong{NUMBER}
        Tolerance, how close to 0 must be.
        NOTE: default value is 10⁻¹⁵}
  \item{n0}{\strong{NUMBER}
        Maximum n0s to run. It will loop n0 times,
        however, if find root, will stop before.}
  \item{VERBOSE}{OPTIONAL \strong{BOOLEAN}
        See results for each iteration}
  \item{VECTOR}{OPTIONAL \strong{BOOLEAN}
        Return a Data frame with calculated values,
        instead of P, i.e., not show if 'Fail'.
  }
  \item{DECIMAL}{OPTIONAL \strong{NUMBER}
        Pass a value to round variables.
        By default, is settled to 13
  }
}
\value{
  \item{DEFAULT}{\strong{FLOAT} Possible solution for a f(x) or 'Fail'.}
  \item{VECTOR=1}{\strong{DATA FRAME} Iterations values, where i is equal to line row}
}
\usage{fixedPointIteration(func, a, b, n0, tol, VERBOSE, VECTOR, DECIMAL)}
\examples{
l <- fixedPointIteration(func='.5 * (10 - x^3)^.5', p0=1.5, tol=10^-9, n0=31, VECTOR=1, DECIMAL=9)
}
\keyword{ root }
\details{Depends of getFuncValue()}
\note{
  \strong{How to see response}
  for (i in 2:nrow(l))

    cat(
      sprintf("\%d\t \%.9f\t \%.9f\t \%.9f\t \%.9f", l[i,1],l[i,2],l[i,3],l[i,4],l[i,5]),
    sep='\n')
}
\seealso{Numerical Analysis 9ed Burden}
\refereces{Numerical Analysis 9ed Burden}
