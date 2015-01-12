\name{TwoStageZ}
\alias{TwoStageZ}
\title{
Zehetmayer's method to analyze data from a two-stage design.
}
\description{
TwoStageZ implements the method proposed by Zehetmayer (2005), which is named as Z method in this package.
}
\usage{
TwoStageZ(x1, s1keep, x2, Vx1, Vx2, gamma1, n1, n2, lambda = 0.5, alpha = 0.05)
}
\arguments{
  \item{x1}{
Observed data from the first stage.
}
  \item{s1keep}{
the list of taxa or gene selected into the second stage.
}
  \item{x2}{
Observed data from the second stage. The order of components, say genes or bacterial species, must be the same as the order used by the first stage data. 
}
  \item{Vx1}{
Variance of first-stage measurements.
}
  \item{Vx2}{
Variance of second-stage measurements.
}
  \item{n1}{
The number of subjects in the first stage.
}
  \item{n2}{
The number of subjects in the second stage.
}
  \item{gamma1}{
P value threshold for the first stage selection of hypotheses into the second stange.
}
  \item{lambda}{
A constant chosen for calculating the null proportion. 0 < lambda < 1. The default is 0.5.
}
  \item{alpha}{
The overall level of FDR. The default is 0.05.
}
}
\value{
  \item{pi0 }{estimated null proportion}
  \item{gamma2 }{estimated second-stage p-value threshold}
  \item{de }{test decisions of each taxa/gene. 1 = significantly differential between the case and control groups; 0 otherwise}
}
\references{
S. Zehetmayer, P Bauer, M. Posch (2005). Two-stage designs for experiments with a large number of hypotheses. Bioinformatics. 21, 3771-3777.
}
\author{
Xiaoshan Wang
}
\seealso{
 \code{\link{TowStageHMM}}
 }
\examples{
NUM1 <- 1000
n1 <- 10
n2 <- 20
z  <- c(rep(0,900),rep(1,100))
Vx1 <- 1
Vx2 <- 1
x1 <- rnorm(NUM1, mean=z,sd=Vx1)

p1 <- 1 - pnorm(as.vector(x1)*sqrt(Vx1),0,1) 
gamma1 <- 0.3
s1keep <- which(p1 < gamma1)

z2 <- z
z2[!s1keep] <- NA

TwoStageZ(x1,s1keep,x2,Vx1,Vx2,gamma1=0.3,n1,n2)
}

