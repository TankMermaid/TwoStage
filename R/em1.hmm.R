em1.hmm <-
function(x, maxiter=200, ptol=tol)
{

NUM<-length(x)
# precision tolerance level
# ptol<-1e-4
niter<-0

### initializing model parameters

pii.new<-c(1, 0)
A.new<-matrix(c(0.8, 0.2, 0.4, 0.6), 2, 2, byrow=T)
f0<-c(0, 1)
f1.new<-c(2, 1)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

pii.old<-pii.new
A.old<-A.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

bwfw.res<-backwardforward1(x, pii.old, A.old, f0, f1.old)

# the hidden states probabilities
gamma<-bwfw.res$pr
# the transition variables
dgamma<-bwfw.res$ts

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

# c. non-null distribution 

q1<-sum(gamma[, 2])
q2<-sum(gamma[, 2]*x)
mu1<-q2/q1
q3<-sum(gamma[, 2]*(x-mu1)*(x-mu1))
sd1<-sqrt(q3/q1)
f1.new<-c(mu1, sd1)

df1<-abs(A.old-A.new)
df2<-abs(f1.old-f1.new)
diff<-max(df1, df2)

}

# f. the final local fdr statistic
lfdr<-gamma[, 1]

# g. return the results of the E-M algorithm

em.var<-list(pii=pii.new, A=A.new, f1=f1.new, lf=lfdr, ni=niter)
return (em.var)

}
