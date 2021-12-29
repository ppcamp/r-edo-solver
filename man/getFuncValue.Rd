\name{getFuncValue}
\alias{getFuncValue}
\title{getFuncValue}
\description{Change x && y variables to actual values. It also evaluates the result.}
\usage{
getFuncValue(func,x,y)
}
\arguments{
  \item{func}{STRING
  Derivative formula write in STRING.
  NOTE: observes that exponential values must use
  '**' instead operator '^', or you must specify
  format=1, to enable author's format.
  }
  \item{x}{Number}
  \item{y}{Number}
}
\details{You must only put right side equation.}
\value{Float value. It returns evaluate values.}
\author{ppcamp}
\examples{getFuncVaue(func="x^2+y-1",x=1,y=2)}
\keyword{ replace }% use one of  RShowDoc("KEYWORDS")
\keyword{ functions }% __ONLY ONE__ keyword per line
