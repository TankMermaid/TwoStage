emC.hmm <-
function(x, L, cons=0.5, maxiter=300, ptol = tol)
{
NUM<-length(x)
# precision tolerance level
# ptol<-1e-3
niter<-0
f0<-c(0, 1)

### initializing model parameters

pii.new<-c(1, 0)
A.new<-matrix(c(0.95, 0.05, 0.5, 0.5), 2, 2, byrow=T)
pc.new<-rep(1, L)/L
mus<-seq(from=-1, by=1.5, length=L)
sds<-rep(1, L)
f1.new<-cbind(mus, sds)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

pii.old<-pii.new
A.old<-A.new
pc.old<-pc.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

bwfw.res<-backwardforward(x, pii.old, A.old, pc.old, f0, f1.old)

# the hidden states probabilities
gamma<-bwfw.res$pr
# the transition variables
dgamma<-bwfw.res$ts
# the weight variables
omega<-bwfw.res$wt

## updating the parameter estimates

# a. initial state distribution

for (i in 0:1)
{
  pii.new[i+1]<-gamma[1, i+1]
}

# b. transition matrix

for (i in 0:1)
{
  for (j in 0:1)
  { 
    q1<-sum(dgamma[i+1, j+1, ])
    q2<-sum(gamma[1:(NUM-1), i+1])
    A.new[i+1, j+1]<-q1/q2  
  }
}

# c. null distribution (given)

# d. non-null mixture distribution

# initializing the vectors of means and variances

mus<-1:L
sds<-1:L

for (c in 1:L)
{

# (i). probability weights
  q1<-sum(omega[, c])
  q2<-sum(gamma[, 2])
  pc.new[c]<-q1/q2
  
# (ii). means
  q3<-sum(omega[, c]*x)
  mus[c]<-q3/q1

# (iii). sds
  q4<-sum(omega[, c]*(x-mus[c])*(x-mus[c]))
  sds[c]<-sqrt(q4/q1)
  if (sds[c] < cons) sds[c] = cons
}
# the non-null mixture distribution
f1.new<-cbind(mus, sds)

diff<-max(abs(f1.old-f1.new))

}

# f. the final local fdr statistic
lfdr<-gamma[, 1]

# g. return the results of the E-M algorithm

em.var<-list(pii=pii.new, A=A.new, pc=pc.new, f1=f1.new, lf=lfdr, ni=niter)
return (em.var)

}
