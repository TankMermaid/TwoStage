\name{TwoStageHMM}
\alias{TwoStageHMM}
\title{
HMM-based FDR control for two-stage design based multiple hypothesis tests under dependence 
}
\description{
The TwoStageHMM implements the hidden Markov model based FDR control as described in Wang et al., 2014 for two-stage design based multiple comparison adjustment under dependence. 
}
\usage{
TwoStageHMM(x1, s1keep, x2, L = 2, maxiter = 100, EM.const = "No", tol = 1e-04, alpha = 0.05)
}
\arguments{
  \item{x1}{
A vector of the first stage observed data.
}
  \item{s1keep}{
the list of taxa or gene selected into the second stage
}
  \item{x2}{
A vector of the second stage observed data.
}

  \item{L}{
the number of components of the normal mixture of f1 in the hidden markov model. The default is 2.
}
  \item{maxiter}{
the maximum number of iteration.
}
  \item{EM.const}{
Choose unconstrained EM algorithm or constrained EM algorithm for HMM model. The default is "NO" for unconstrained EM algorithm.
}
  \item{tol}{
Tolerance for convergence. The default is 1e-10.
}
  \item{alpha}{
the nominal level of FDR. The Default is 0.05.
}
}
\value{
  \item{nr }{the total number of rejections among the retested hypotheses.}
  \item{re }{the list of rejected taxa/gene}
  \item{ac }{the list of accepted taxa/gene}
  \item{de }{test decisions of each taxa/gene. 1 = significantly differential between the case and control groups; 0 otherwise.}
}
\references{
X Wang, J. Starr, F. Jorge (2014) Two-stage design optimizes power and false discovery rates for microbiome research.

W. Sun and T. Cai (2009). Large-scale multiple testing under dependence. J. R.Stat. Soc B. 71, 393-424.
}
\author{
Xiaoshan Wang
}
\seealso{
 \code{\link{TwoStageZ}}
}
\examples{
################
# Simulate data
################
NUM1 = 1000
pii<-c(1, 0)  
f0<-c(0, 1)
f1<-c(2, 1)
f2 <-c(2.5,1)

a11 <- 0.8
A<-matrix(c(0.95, 0.05, 1- a11, a11), 2, 2, byrow=T)

dt1<-sim1(NUM1, pii, A, f0, f1, f2)

x1 <- dt1$o1
x2 <- dt1$o2

theta <- dt1$s
s1kp <- which(theta==1)
xs2 <- x2[s1kp]

two <- TwoStageHMM(x1, s1kp, xs2, L=2, maxiter=500, EM.const = 'No')
}
