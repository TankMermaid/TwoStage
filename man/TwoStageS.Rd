\name{TwoStageS}
\alias{TwoStageS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simes' combination function based analysis of two-stage data.
}
\description{
TwoStageS implements the Simes' combination method proposed by Sarkar (2013), which combines two-stage p values by Simes' combination function.
}
\usage{
TwoStageF(p1, p2, alpha = 0.05, alpha.L0 = 0, alpha.U0 = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
Stage 1 p values.
}
  \item{p2}{
Stage 1 p values. Require to have the exact same order of Stage 1 p values, with NA replacing excluded tests by Stage 1 screening.
}
  \item{alpha}{
Overall FDR level.
}
  \item{alpha.L0}{
A FDR level set for declaring significant test in Stage 1. Default is 0. Its value must be smaller than overall FDR alpha.
}
  \item{alpha.U0}{
Stage 1 threshold to screen tests. p values between alpha.L0 and alpha.U0 will be re-tested in Stage 2.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Number of rejection and acceptance }{Includes overall FDR alpha, Stage 1 FDR spending of alpha.L0, Stage 1 screening threshold alpha.U0, Stage 1 number of rejection R1, the number to be retested S1, Stage 1 number of acceptance A1,Stage 1 number of rejection R2, overall number of rejection R.}
  \item{List of rejected test }{A list of the order of rejected tests.}
%% ...
}
\references{
SK Sarkar, J. Chen, W. Guo (2013). Multiple Testing in a Two-Stage Adaptive Design with Combination Tests Controlling FDR. J. Am. Stat. Assoc. 108: 1385-1401.
}
\author{
Adapted from the codes of J. Chen, the author of the reference paper. 
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{TowStageHMM}},\code{\link{TowStageZ}}
}
\examples{
p1 <- runif(100)
p2 <- ifelse(p1 < 0.2,runif(1,0,0.1),NA)
TwoStageF(p1, p2, alpha = 0.05, alpha.L0 = 0, alpha.U0 = 0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
