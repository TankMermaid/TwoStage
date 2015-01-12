sim <-
function(NUM, pii, A, f0, pc, f1, f2)
{

theta<-rep(0, NUM)
x1 <- rep(0, NUM)
x2 <- rep(0, NUM)
nc<-length(pc)

if (nc==1)
{ 
  data<-rdata1.hmm (NUM, pii, A, f0, f1, f2)
}    

else
{
theta[1]<-rbinom(1, 1, pii[2])
for (i in 2:NUM)
{
  if (theta[i-1]==0)
     theta[i]<-rbinom(1, 1, A[1, 2])
  else
     theta[i]<-rbinom(1, 1, A[2, 2])
}

## generating stage-1 observations
for (i in 1:NUM)
{
  if (theta[i]==0)
  {
    mu0<-f0[1]
    sd0<-f0[2]
    x1[i]<-rnorm(1, mean=mu0, sd=sd0)
  }
  else
  { 
    c<-sample(1:nc, 1, prob=pc)
    mu1<-f1[c, 1]
    sd1<-f1[c, 2]
    x1[i]<-rnorm(1, mean=mu1, sd=sd1)
  }
}

## generating stage-2 observations
for (i in 1:NUM)
{
  if (theta[i]==0)
  {
    mu0<-f0[1]
    sd0<-f0[2]
    x2[i]<-rnorm(1, mean=mu0, sd=sd0)
  }
  else
  { 
    c<-sample(1:nc, 1, prob=pc)
    mu1<-f2[c, 1]
    sd1<-f2[c, 2]
    x2[i]<-rnorm(1, mean=mu1, sd=sd1)
  }
}

data<-list(s=theta, o1=x1, o2=x2)
}

return (data)

}
