TwoStageHMM <-
function(x1, s1keep, x2, L=2, maxiter=100, EM.const = 'No', tol=1e-4, alpha=0.05){

if(EM.const == 'No'){
 cat('Running with HMM, Normal components',L,', Unconstrained EM','...','\n')
 if (L==2)  pc = 1
 res1 <- try(em1.hmm(x=x1, ptol = tol, maxiter=maxiter))
 res2 <- try(em1.hmm(x=x2, ptol = tol, maxiter=maxiter))
 lsi = cbind(res1$lf[s1keep],res2$lf)
 lsi.m = apply(lsi,1,mean)
 res<-mt.hmm(lsi.m, alpha)
}

if(EM.const != 'No'){
 cat('Running with HMM, Normal components',L,', Constrained EM','...','\n')
 if (L==2)  pc = 1
 res1 <- try(em1C.hmm(x=x1, ptol = tol, maxiter=maxiter))
 res2 <- try(em1C.hmm(x=x2, ptol = tol, maxiter=maxiter))
 lsi = cbind(res1$lf[s1keep],res2$lf)
 lsi.m = apply(lsi,1,mean)
 res<-mt.hmm(lsi.m, alpha)
}

return(res)
}
