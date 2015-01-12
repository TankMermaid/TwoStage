powz <-
function(gamma1, r=0.5,m1=1000, pi0=0.99,alpha=0.05,delta=1){

N = 8*m1
n1 = N*r/m1

pow1 = 1 - pnorm(qnorm(gamma1)-delta*sqrt(n1/2))

m2 = m1*(pi0*gamma1 + (1-pi0)*pow1)
n2 = N*(1-r)/m2
nt = n1 + n2

gfun<-function(z,g2) (1-pnorm((qnorm(1-g2)-sqrt(n1/nt)*z)/sqrt(n2/nt)))*dnorm(z)

gam<-function(g2) integrate(gfun,lower=qnorm(1-gamma1),upper=Inf,g2=g2) 

pfun <- function(z,g2) (1-pnorm((qnorm(1-g2)-sqrt(n1/nt)*z)/sqrt(n2/nt),sqrt(n2)*delta,1))*dnorm(z,sqrt(n1)*delta,1)

pow <- function(g2) integrate(pfun,lower=qnorm(1-gamma1),upper=Inf,g2=g2)

fdr <- function(g2) pi0*gam(g2)$value/(pi0*gam(g2)$value + (1-pi0)*pow(g2)$value)

f.fdr = function(g2) fdr(g2) - alpha

g2h = uniroot(f.fdr,c(0.00000001,0.01),tol = 1e-200)$root

opow = pow(g2h)$value

return(opow)
}
