\name{plotPatterns}
\alias{plotPatterns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot methylation pattern distribution
}
\description{
  Plot the observed distribution and the estimated distribution of the methylation patterns
}
\usage{
  plotPatterns(compareData, yLimit1=NULL, yLimit2=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compareData}{
    data frame, obtained from the output of the function \code{\link{estimatePatterns}}
  }
  \item{yLimit1}{
    upper limit of y-axis on left hand scale of the first graph. If NULL, defaults to show all patterns
  }
  \item{yLimit2}{
    upper limit of y-axis on left hand scale of the second graph. If NULL, defaults to show most patterns
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  estimates <- estimatePatterns(patternsExample, 
                                epsilon=0.02, 
                                eta=0.01)
  plotPatterns(estimates[[1]])
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

% vim:ts=2:sw=2:sts=2:expandtab:
