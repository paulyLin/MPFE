\name{patternMap}
\alias{patternMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot a representation of the patterns and their frequencies
}
\description{
  Plot the observed distribution and the estimated distribution of the methylation patterns
}
\usage{
  patternMap(patterns,
             minFreq=0,
             maxFreq=1,
             noSpurious=TRUE,
             estimatedDistribution=TRUE,
             topDown=TRUE,
             allTicks=FALSE,
             methCol='black',
             unMethCol='white',
             ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patterns}{
    A data frame obtained from the output of the function \code{\link{estimatePatterns}}.
  }
  \item{minFreq}{
    Only plot patterns with at least \code{minFreq} frequency.
  }
  \item{maxFreq}{
    Only plot patterns with more \code{maxFreq} frequency or more.
  }
  \item{noSpurious}{
    Don't plot spurious patterns (only relevant if \code{estimatedDistribution} is \code{FALSE}).
  }
  \item{estimatedDistribution}{
    Use the frequencies from the estimated distribution. If \code{FALSE}, use the observed distribution.
  }
  \item{topDown}{
    Put the most abundant patterns at the top. If \code{FALSE} the most abundant patterns are at the bottom.
  }
  \item{allTicks}{
    Draw a tick under every position.
  }
  \item{methCol}{
    The colour for the methylated positions. Can be a single colour, a vector of colours (recycled), or a function (for instance from \code{\link{colorRampPalette}}).
  }
  \item{unMethCol}{
    As \code{methCol} but for un-methylated positions.
  }
  \item{...}{
    Other  arguments passed to \code{\link{plot}}
  }
}
\details{
	This function draws a map of the different pattern and their frequencies based on the values returned by \code{\link{estimatePatterns}}.
}

%% \value{
%%	A plot that compares the observed read distribution with the estimated distribution. 
%% }

\author{
  Peijie Lin, Sylvain Foret, Conrad Burden
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
  data(patternsExample)
  estimates <- estimatePatterns(patternsExample, 
                                epsilon=0.02, 
                                eta=0.01)
  patternMap(estimates[[1]])
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

% vim:ts=2:sw=2:sts=2:expandtab:
