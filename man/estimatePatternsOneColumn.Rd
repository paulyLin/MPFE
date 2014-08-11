\name{estimatePatternsOneColumn}
\alias{estimatePatternsOneColumn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate distribution of methylation patterns with one column of counts
}
\description{
  This function is part of the inner working of the package. 
Please refer to \code{\link{estimatePatterns}}
}
\usage{
  estimatePatterns(patternCounts,
                   epsilon,
                   eta,
                   column,
                   fast,
                   steps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patternCounts}{
    data frame with methylation patterns in first column and pattern counts in
    subsequent columns.
  }
  \item{epsilon}{
    non-converson rate, a value between 0 and 1.
  }
  \item{eta}{
    error rate, either a vector of numbers between 0 and 1 of length equal to
    the number of CpG sites or a single value between 0 and 1 for a single error
    rate across all sites.
  }
  \item{column}{
    index of the column of `patternCounts' to process, an integer value from 1
    to the number of pattern counts columns in `patternCounts'.  
  }
  \item{fast}{
    logical, if TRUE, fast version implemented (default).
  }
  \item{steps}{
    number of steps for the optimiser, passed to \code{\link{constrOptim}}.
    If NULL, defaults to 20000 steps.
  }
}

\value{
The function returns a data frame,
which contains the following columns:
  \item{pattern}{
    the list of input patterns (factor)
  }
  \item{coverage}{
    the number of reads for each pattern (integer)
  }
  \item{observedDistribution}{
    the observed frequencies of each pattern (numeric)
  }
  \item{estimatedDistribution}{
    the estimated frequencies (numeric)
  }
  \item{spurious}{
    indicates whether the patterns are real or spurious (logical)
  }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Peijie Lin, Sylvain Foret, Conrad Burden
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  data(patternsExample)
  estimatePatternsOneColumn(patternsExample, 
                   	    epsilon=0.02, 
                   	    eta=0.01,
                   	    column=1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

% vim:ts=2:sw=2:sts=2:expandtab:
