\name{odeEuler}
\alias{odeEuler}
\title{odeEuler}
\description{Calculates ordinary differentials equations using Euler's method}
\usage{odeEuler(func,a,b,m,y0,VERBOSE)}
\arguments{
  \item{func}{\strong{STRING}
        Derivative formula write in STRING.}
  \item{a}{\strong{INT}
        Lower boundary}
  \item{b}{\strong{INT}
        Upper boundary}
  \item{m}{\strong{INT}
        Amount of intervals}
  \item{y0}{\strong{INT}
        Start value}
  \item{VERBOSE}{\emph{OPTIONAL} \strong{BOOLEAN}
        See results for each iteration}
  \item{VECTOR}{OPTIONAL \strong{BOOLEAN}
        Return a Data frame with calculated values,
        instead of P, i.e., not show if 'Fail'.}
}
\value{It returns a \strong{VECTOR c()} with x,y, && fxy values referents to m iterations.
  \item{x}{ \strong{FLOAT}}
  \item{y}{ \strong{FLOAT}}
  \item{fxy}{ \strong{FLOAT}}
  Or then if
  \item{VECTOR=1}{\strong{DATA FRAME} Iterations values, where i is equal to line row}
}
\author{@ppcamp}
\note{
Pseudocode (p. 295)

Author: Frederico Ferreira Filho

Title: Algoritmos Numéricos 2ª}

\examples{
odeEuler(func='x-2*y+1',a=0,b=1,y0=1,m=10, VERBOSE=1)
odeEuler(func='y-x^2+1',a=0,b=2,y0=0.5,m=10,VERBOSE=1)}
\keyword{ derivative }% use one of  RShowDoc("KEYWORDS")
