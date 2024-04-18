#'@title DOTNB_pdf
#'@description
#'Calculates the p-value for the probability of significant difference between two negative binomial distributions using the DOTNB probability density function.
#'
#'@param r1 The r parameter for the first distribution
#'@param p1 The p parameter for the first distribution
#'@param r2 The r parameter for the second distribution
#'@param p2 The p parameter for the second distribution
#'@param tt the observed difference between the two distributions
#'@return the PDF output
#'@export

DOTNB_pdf<-function(r1,p1,r2,p2,tt){
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

#'@title DOTNB_mean
#'@description
#'Calculates the mean for the DOTNB distribution based on the parameters of two negative binomial distributions
#'
#'@param r1 The r parameter for the first distribution
#'@param p1 The p parameter for the first distribution
#'@param r2 The r parameter for the second distribution
#'@param p2 The p parameter for the second distribution
#'@return The mean for the representative DOTNB distribution
#'@export
DOTNB_mean<-function(r1,p1,r2,p2){
  q1 = 1-p1
  q2 = 1-p2
  dmean = ((r1*q1)/p1) - (r2*q2/p2)
  return(dmean)
}

#'@title DOTNB_mean
#'@description
#'Calculates the variance for the DOTNB distribution based on the parameters of two negative binomial distributions
#'
#'@param r1 The r parameter for the first distribution
#'@param p1 The p parameter for the first distribution
#'@param r2 The r parameter for the second distribution
#'@param p2 The p parameter for the second distribution
#'@return The variance for the representative DOTNB distribution
#'@export
DOTNB_var<-function(r1,p1,r2,p2){
  q1 = 1-p1
  q2 = 1-p2
  dvar = (r1*q1)/(p1^2) + (r2*q2)/(p2^2)
  return(dvar)
}

#'@title DOTNB_pdf
#'@description
#'Calculates the p-value for the probability of significant difference between two negative binomial distributions using the DOTNB cumulative density function.
#'
#'@param r1 The r parameter for the first distribution
#'@param p1 The p parameter for the first distribution
#'@param r2 The r parameter for the second distribution
#'@param p2 The p parameter for the second distribution
#'@param tt the observed difference between the two distributions
#'@return the CDF output
#'@export
DOTNB_cdf<-function(r1,p1,r2,p2,tt){
  cvalue=0
  pdfv=DOTNB_pdf(r1,p1,r2,p2,tt)
  maxiter = 2000
  if(is.infinite(pdfv)){
    if (tt > DOTNB_mean(r1,p1,r2,p2)){
      cvalue=1
    }
  }else{
    iter = 1
    while(is.finite(pdfv)){
      cvalue=cvalue+pdfv
      tt=tt-1
      pdfv=DOTNB_pdf(r1,p1,r2,p2,tt)
      iter = iter +1
      if(iter > maxiter){
        break
      }
    }
  }
  return(cvalue)
}

