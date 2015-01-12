sim1 <-
function(NUM, pii, A, f0, f1, f2)
{

theta<-rep(0, NUM)
x1 <- rep(0, NUM)
x2 <- rep(0, NUM)

## generating the states
 # initial state
theta[1]<-rbinom(1, 1, pii[2])
 # other states
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
    x1[i]<-rnorm(1, mean=f0[1], sd=f0[2])
  }
  else
  { 
    x1[i]<-rnorm(1, mean=f1[1], sd=f1[2])
  }
}

## generating stage-2 observations

for (i in 1:NUM)
{
  if (theta[i]==0)
  {
    x2[i]<-rnorm(1, mean=f0[1], sd=f0[2])
  }
  else
  { 
    x2[i]<-rnorm(1, mean=f2[1], sd=f2[2])
  }
}

data<-list(s=theta, o1=x1, o2=x2)
return (data)
}
