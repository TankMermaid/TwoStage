TwoStageZ <-
function(x1,s1keep,x2,gamma1,n1,n2,lambda=0.5,alpha=0.05){

xx2 = rep(NA, length(x1))
xx2[s1keep] = x2
if (length(Vx1)==1) v1 = rep(Vx1, length(x1))
else v1 = Vx1
if (length(Vx2)==1) v2 = rep(Vx2, length(x1))
else v2 = Vx2

z = x1
p1 <- 1 - pnorm(as.vector(x1)*sqrt(v1),0,1)
p <- 1 - pnorm((x1*v1/sqrt(v1+v2) + xx2*v2/sqrt(v1+v2)),0,1)

pi0.est<-function(p1,lambda) min(sum(p1>lambda)/((1-lambda)*length(p1)),1)

gfun<-function(z,n1,n2,g2) (1-pnorm((qnorm(1-g2)-sqrt(n1/(n1+n2))*z)/sqrt(n2/(n1+n2))))*dnorm(z)

int<-function(gamma1,n1,n2,g2) integrate(gfun,lower=qnorm(1-gamma1),upper=Inf,n1=n1,n2=n2,g2=g2) 

pi0 = pi0.est(p1,lambda)

fdr<-function(p,p1,gamma1,n1,n2,g2,lambda) pi0*length(p1)*int(gamma1,n1,n2,g2)$value/(max(sum(p<g2,na.rm=TRUE),1))
    
fdrminus<-function(p,p1,gamma1,n1,n2,g2,lambda,alpha) fdr(p,p1,gamma1,n1,n2,g2,lambda)-alpha

gamma2<-function(p,p1,gamma1,n1,n2,lambda,alpha) uniroot(fdrminus,c(0.00000001,0.05),tol = 1e-200,p=p,p1=p1,gamma1=gamma1,n1=n1,n2=n2,lambda=lambda,alpha=alpha)

gam2 = gamma2(p,p1,gamma1,n1,n2,lambda,alpha)$root

testdecision<-function(p,p1,gamma1,n1,n2,lambda,alpha) (p1<=gamma1) & (p<gam2)
    
tdec = testdecision(p,p1,gamma1,n1,n2,lambda,alpha)

ones = which(tdec==T)
zeros = rep(0,length(p1))
zeros[ones] = 1
out = list(pi0 = pi0, gamma2 = gam2, de=zeros)
return(out)
}
