#!/usr/bin/env python3

from scipy.special import gamma, factorial, hyp2f1

"""
This script provides DOTNB functions
Written by Ruoyu Chen
"""

def DOTNB_mean(r1, p1, r2, p2, k):
    q1 = 1-p1
    q2 = 1-p2
    dmean = ((r1*q1)/p1) - (r2*q2/p2)
    return dmean


def DOTNB_var(r1, p1, r2, p2, k):
    q1 = 1-p1
    q2 = 1-p2
    dvar = ((r1*q1)/p1**2) - (r2*q2/p2**2)
    return dvar

def DOTNB_pdf(r1, p1, r2, p2, k):
    q1 = 1-p1
    q2 = 1-p2
    d = 0
    if k > 0:
        d = (p1**r1)*(p2**r2)*((gamma(r1+k)/gamma(r1))/gamma(1+k))*(q1**k)*hyp2f1(r1+k,r2,k+1,q1*q2)
    else:
        nk = -k
        d = (p1**r1)*(p2**r2)*((gamma(r2+k)/gamma(r2))/gamma(1+nk))*(q2**nk)*hyp2f1(r2+nk,r1,nk+1,q1*q2)
    return d

def DOTNB_cdf(r1, p1, r2, p2, k):
    cvalue = 0
    pdfv = DOTNB_pdf(r1, p1, r2, p2, k)
    TOLERANCE = 1.0e-20
    i = 0
    maxiter = 2000
    
    if pdfv < TOLERANCE:
        if pdfv > DOTNB_mean(r1, p1, r2, p2):
            cvalue = 1
        else:
            cvalue = 0
    else:
        while (pdfv > TOLERANCE):
            cvalue = cvalue + pdfv
            k = k-1
            i = i + 1
            pdfv = DOTNB_pdf(r1, p1, r2, p2, k)
            if i > maxiter:
                break
                
    return cvalue
