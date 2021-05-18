diff_pvalue<-function(r1,p1,r2,p2){
  pvalue=min(1-DOTNB_cdf(r1,p1,r2,p2,0),DOTNB_cdf(r1,p1,r2,p2,-1))
  return(pvalue)
}