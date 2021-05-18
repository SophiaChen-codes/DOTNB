
DOTNB_pdf<-function(r1,p1,r2,p2,tt){
# calculate the p-value of two different NB distributions (r1,p1) and (r2,p2)
# tt is the observed number of difference
  
  if(!require('hypergeo',character.only = TRUE)) install.packages('hypergeo')
  library('hypergeo',character.only = TRUE)
  
  q1 = 1-p1
  q2 = 1-p2
  if(tt>0){
    dd = (p1^r1)*(p2^r2)*((gamma(r1+tt)/gamma(r1))/gamma(1+tt))*(q1^tt)*hypergeo(r1+tt,r2,tt+1,q1*q2);
  }else{
    ntt=-tt
    dd = (p1^r1)*(p2^r2)*((gamma(r2+ntt)/gamma(r2))/gamma(1+ntt))*(q2^ntt)*hypergeo(r2+ntt,r1,ntt+1,q1*q2);
  }
  return(Re(dd))
}

DOTNB_pdf2<-function(r1,p1,r2,p2,tt){
# calculate the p-value of two different NB distributions (r1,p1) and (r2,p2)
# tt is the observed number of difference
  
  if(!require('hypergeo',character.only = TRUE)) install.packages('hypergeo')
  library('hypergeo',character.only = TRUE)
  
  q1 = 1-p1
  q2 = 1-p2
  if(tt>0){
    dd = (p1^r1)*(p2^r2)*((gamma(r1+tt)/gamma(r1))/gamma(1+tt))*(q1^tt)*genhypergeo(c(r1+tt,r2),tt+1,q1*q2);
  }else{
    ntt=-tt
    dd = (p1^r1)*(p2^r2)*((gamma(r2+ntt)/gamma(r2))/gamma(1+ntt))*(q2^ntt)*genhypergeo(c(r2+ntt,r1),ntt+1,q1*q2);
  }
  return(dd)
}

DOTNB_mean<-function(r1,p1,r2,p2){
# calculate the mean of the difference of two NB distributions (r1,p1) and (r2,p2)
  
  q1 = 1-p1
  q2 = 1-p2
  dmean = ((r1*q1)/p1) - (r2*q2/p2)
  return(dmean)
}

DOTNB_var<-function(r1,p1,r2,p2){
# calculate the variance of the difference of two NB distributions (r1,p1) and (r2,p2)
  q1 = 1-p1
  q2 = 1-p2
  dvar = (r1*q1)/(p1^2) + (r2*q2)/(p2^2)
  return(dvar)
}

DOTNB_cdf<-function(r1,p1,r2,p2,tt){
# calculate the cumulative probability of tt that is the difference of two NB distributions (r1,p1) and (r2,p2)
# tt is the observed number of difference
# sum pvalue for <= tt, 
  cvalue=0
  pdfv=DOTNB_pdf(r1,p1,r2,p2,tt)
  if(is.infinite(pdfv)){
    if (tt > DOTNB_mean(r1,p1,r2,p2)){
      cvalue=1
    }
  }else{
    while(is.finite(pdfv)){
      cvalue=cvalue+pdfv
      tt=tt-1
      pdfv=DOTNB_pdf(r1,p1,r2,p2,tt)
    }
  }
  return(cvalue)
}

