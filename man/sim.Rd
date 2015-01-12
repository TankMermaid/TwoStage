\name{sim}
\alias{sim}
\title{
Data simulation for hidden Markov model.
}
\description{
Generate data based on predefined transition probability matrix.
}
\usage{
sim(m, pii, A, f0, pc, f1, f2)
}
\arguments{
  \item{m}{
The totmal number of hypothesis.
}
  \item{pii}{
Initial status. The default is pii = (0,1).
}
  \item{A}{
The transition probability matrix.
}
  \item{f0}{
The distribution to generate data of null hypotheses. 
}
  \item{pc}{
The proportions for non-null components. For example, if we simulate 0.3N(1,1) + 0.3N(2,1) + 0.4N(3,1), then pc = (0.3, 0.3, 0.4).
}
  \item{f1}{
The distribution to generate data of non-null hypotheses for the first stage.  It should be a vector for one-component non-null, and a matrix with means of non-null components on Column 1 and corresponding standard deviations on Colum 2.
}
  \item{f2}{
The distribution to generate data of non-null hypotheses for the second stage. It's struction should be the same as f1.
}
}
\value{
 \item{s }{The true status.}
 \item{o1 }{Generated data for the first stage.}
 \item{o2 }{Generated data for the second stage.}
}

\examples{
m = 1000
pii<-c(1, 0)  
f0<-c(0, 1)
f1<-matrix(c(2, 1, 2.5, 1),ncol=2, byrow=T)
f2<-matrix(c(2, 1, 2.5, 1),ncol=2, byrow=T)
a11 = 0.8
A<-matrix(c(0.95, 0.05, 1- a11, a11), 2, 2, byrow=T)
pc = c(0.5, 0.5)

dt<-sim(m, pii, A, f0, pc, f1, f2)
}


