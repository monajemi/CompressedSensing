\name{PT-package}
\alias{PT-package}
\alias{PT}
\docType{package}
\title{
A package for analysis of Phase Transition data
}
\description{
Various functions for analysis of compressed
sensing phase transition data
}
\details{
\tabular{ll}{
Package: \tab PT\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-01-27\cr
License: \tab LGPL >= 2.0 \cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Hatef Monajemi
Maintainer: Hatef Monajemi (monajemi@stanford.edu)
}
\references{
Deterministic matrices matching the compressed sensing phase transitions of Gaussian random matrices.
http://www.pnas.org/content/110/4/1181.short
}

\keyword{ PT }
\seealso{

}
\examples{
delta = 1/2;
field = "R";
type  = "w";
epsilon <- predictPT(delta, field, type);
}
