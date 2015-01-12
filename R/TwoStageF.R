TwoStageF <-
function(p1,p2,alpha=0.05,alpha.L0=0,alpha.U0=0.3){


matchfirst <- function(a, v) {match(a, v, nomatch=length(v)+1)}


m<-length(p1)   ## total number of hypothese to test
n<-1       ## just 1 iteration


## stage 1 p-value via t-test
## order p1 for stage 1 stepup and stepdown tests

p1<-matrix(p1,nrow=n,ncol=m)
p1.up<-t(apply(p1,1,sort,decreasing=FALSE))
p1.down<-t(apply(p1,1,sort,decreasing=TRUE))


## find early rejection and early acceptance
#### find number of early rejections and early acceptance at stage 1 
 

   FDR1<-c(rep(0,n*m))
   FDR1<-matrix(FDR1,nrow=n,ncol=m)

   reject1<-c(rep(0,n*m))
   reject1<-matrix(reject1,nrow=n,ncol=m)


## calculate FDR1 using alternative definition


for (i in 1:n){
   for (k in 1:m) {
    FDR1[i,k]=(m*p1.up[i,k])/k
   }
}



   ## find FDR1<=alpha.L0 using stepdown procedure
   ## FDR1 vs alpha.L0 find the first FDR1(i)>alpha.L0 then accept H(i)~H(m)

   compare.L1<-c(rep(0,m*n))
   compare.L1<-matrix(compare.L1,nrow=n,ncol=m)

   compare.L1<-(FDR1>alpha.L0)


   ## find the first acceptance where a=0, 

   R1<-matrix(0,nrow=n,ncol=1)

   accept<-apply(compare.L1, 1, matchfirst, a=1)
   R1<-accept-1


  ## find FDR1>alphal.U0 using stepup procedure, FDR1.down
  ## FDR1 vs alpha.U0 find the first FDR1(i)<alpha.U0 then accept H(i)~H(m)

   compare.U1<-c(rep(0,m*n))
   compare.U1<-matrix(compare.U1,nrow=n,ncol=m)   
   
   
   FDR1.new<-matrix(0,nrow=n,ncol=m)  

for (i in 1:n){
 for (j in 1:m) {
      FDR1.new[i,j]<-FDR1[i,(m-j+1)]
   }
  }

   compare.U1<-(FDR1.new<=alpha.U0)

   reject<-apply(compare.U1, 1, matchfirst, a=1)

   A1<-reject-1
   S1<-m-A1


#print(R1)
#print(S1)
 



## stage 2 p-value via t-test, 


p2<-matrix(p2,nrow=n,ncol=m)


## get combined p-value via fisher's combination function

m2<-S1-R1


p<-c(rep(0,m*n))   ## combined p-value
p<-matrix(p1*p2,nrow=n,ncol=m)
p[is.na(p)] = 2

#for (i in 1:n){
#  for (k in 1:m){
# if (p1[i,k]>(alpha.L0*R1[i]/m) & p1[i,k]<=(alpha.U0*S1[i]/m)) {p[i,k]<-p1[i,k]*p2[i,k]} else {p[i,k]<-2}
# }
#}

p.order<-t(apply(p,1,sort,decreasing=FALSE))
p.ord = order(p,decreasing = F)


## apply two-stage BH-TSADC procedure with fisher's combination function, with pi0hat estimate

###########################################
### control FDR2 at alpha-alpha.L0    #####
###########################################


 pi0.hat<-(m-S1+1)/(m*(1-alpha.U0))


   FDR2<-c(rep(0,n*m))
   FDR2<-matrix(FDR2,nrow=n,ncol=m)

   reject2<-c(rep(0,n*m))
   reject2<-matrix(reject2,nrow=n,ncol=m)


## calculate FDR2 using alternative definition


for (i in 1:n){
   for (k in 1:m) {
      if (k<=S1[i]-R1[i]){

if (p.order[i,k]<(R1[i]*alpha.L0)/m) {FDR2[i,k]=(m*pi0.hat[i]*p.order[i,k]*log((S1[i]*alpha.U0)/(R1[i]*alpha.L0)))/(R1[i]+k)}
if (p.order[i,k]>=(R1[i]*alpha.L0)/m & p.order[i,k]<(S1[i]*alpha.U0)/m) {FDR2[i,k]=m*pi0.hat[i]*(p.order[i,k]-(R1[i]*alpha.L0)/m+p.order[i,k]*log((S1[i]*alpha.U0)/(m*p.order[i,k])))/(R1[i]+k)}
if (p.order[i,k]>=(S1[i]*alpha.U0)/m) {FDR2[i,k]=pi0.hat[i]*(S1[i]*alpha.U0-R1[i]*alpha.L0)/(R1[i]+k)}

        }
      else {FDR2[i,k]<-2}
   }
}


  ## find FDR2>alpha-alpha.L0 using stepup procedure, FDR2.down
  ## FDR2 vs alpha-alpha.L0 find the first FDR2(i)<alpha-alpha.L0 then accept H(i)~H(m)

## get the number of rejections at stage 2 

     compare.c<-c(rep(0,m*n))
     compare.c<-matrix(compare.c,nrow=n,ncol=m)   
   

FDR2.new<-matrix(0,nrow=n,ncol=m)  

for (i in 1:n){
 for (j in 1:m) {

      FDR2.new[i,j]<-FDR2[i,(m-j+1)]

   }
  }
   
   compare.c<-(FDR2.new<=(alpha-alpha.L0))

   reject2<-apply(compare.c, 1, matchfirst, a=1)

   A2<-reject2-1
   R2<-m-A2
  
  #print(R2)


R<-R1+R2

compare.c2 = as.vector((FDR2<=(alpha-alpha.L0))*1)
pos.list = p.ord[which(compare.c2==1)]

result<-data.frame(alpha, alpha.L0, alpha.U0,R1,S1,A1,R2,R)
colnames(result)<-c("alpha", "alpha.L0", "alpha.U0","R1","S1","A1","R2","R")

return(list(result,pos.list))
}
