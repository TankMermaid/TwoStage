PowerZ <-
function(m, pi0, delta){
 rs=seq(0.01,1,by=0.01)
 pows = numeric()
 for (j in 1:length(rs)){
  pows[j] = optimize(powz,c(0,0.5),maximum =T, r = rs[j], m1=m, pi0 = pi0, delta=delta)$objective
 }

 #if ("plot") plot(rs, pows)
 rr = rs[which(pows==max(pows))]
 ropt = optimize(powz,c(0,0.5),maximum =T, r = rr)
 pow.rr = ropt$objective
 g1.rr = ropt$maximum 
 return(c(gamma1 = g1.rr, n.ratio=rr, max.power=pow.rr))
}
