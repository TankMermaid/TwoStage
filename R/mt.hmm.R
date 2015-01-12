mt.hmm <-
function(lfdr, q)
{
  m=length(lfdr)
  st.lfdr<-sort(lfdr)
  o.lfdr = order(lfdr,decreasing=F)
  
  hps<-rep(0, m)
  if (min(lfdr)>q)
  {
    k<-0
    threshold<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k=1
    while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
      k=k+1
    }
    k<-k-1
    
    # threshold<-st.lfdr[k]
    reject<-o.lfdr[1:k] 
    accept<-o.lfdr[(k+1):m]
    hps[reject]<-1
  }
  y<-list(nr=k, re=reject, ac=accept, de=hps)
  return (y)

}
