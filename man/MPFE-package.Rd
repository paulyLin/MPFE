\name{MPFE-package}
\alias{MPFE-package}
\alias{MPFE}
\docType{package}
\title{
  MPFE
}

\description{
  Estimate distribution of methylation patterns from a table of counts from a
  bisulphite sequencing experiment given a non-conversion rate and sequencing
  error rate.
}

\details{
  \tabular{ll}{
    Package: \tab MPFE\cr
    Type: \tab Package\cr
    License: \tab GPL(>=3)\cr
  }
  The main component of this package is the function \code{\link{estimatePatterns}}, 
  which reads a table of read counts of bisulphite sequencing data for a given
  amplicon and generates a table and plot of the estimated distribution over
  methylation patterns.
}

\author{
  Peijie Lin, Sylvain Foret, Conrad Burden

  Maintainer: conrad.burden@anu.edu.au
}


\keyword{
  amplicon,
  bisulphite sequencing,
  bisufite sequencing,
  methylation
}


\examples{
  data(patternsExample)
  estimates <- estimatePatterns(patternsExample, epsilon=0.02, eta=0.01)
  estimates
  plotPatterns(estimates[[2]])
}

% vim:ts=2:sw=2:sts=2:expandtab:
