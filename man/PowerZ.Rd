\name{PowerZ}
%\alias{powerz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size estimation for Zehetmayer's approach.
}
\description{
PowerZ implement a search for optimal allocation of sample size in a two stage design planned to apply Zehetmayer's analysis approach.
}
\usage{
PowerZ(m, pi0, alpha = 0.05, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{gamma1}{
%   Stage 1 p-value threshold to screen hypotheses. Any hypothesis with a p value less than gamma1 will be retested in Stage two.
%}
  \item{m}{
 The total number of hypotheses to be tested.
 }
  \item{r}{
the ratio of sample size allocation between Stage 1 and Stage 2.
}
  \item{alpha}{
 The overall alpha level. Default is 0.05.
}
  \item{pi0}{
 The proportion of null hypotheses.
 }

  \item{delta}{
The effect size between two groups.
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A vector includes Stage 1 threshold gamma1, total sample size n, Stage 1 proportion of sample size r, and the maximum power max.power.
}
\references{
S. Zehetmayer, P Bauer, M. Posch (2005). Two-stage designs for experiments with a large number of hypotheses. Bioinformatics. 21, 3771-3777.
}
\author{
Xiaoshan Wang
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
PowerZ(m=1000,alpha=0.05,pi0=0.95,delta=1.5)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
