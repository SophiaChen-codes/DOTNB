#!/usr/bin/env python3

import scipy.special import gamma, factorial, hyp2f1

"""
This script provides DTONB functions of PDF, mean and variance
"""


def DOTNB_pdf(r1, p1, r2, p2, k):
    q1 = 1-p1
    q2 = 1-p2
    d = 0
    if if k > 0:
        d = (p1**r1)*(p2**r2)*((gamma(r1+k)/gamma(r1))/gamma(1+k))*(q1**k)*hyp2f1(r1+k,r2,k+1,q1*q2)
    else:
        nk = -k
        d = (p1**r1)*(p2**r2)*((gamma(r2+k)/gamma(r2))/gamma(1+nk))*(q2**nk)*hyp2f1(r2+nk,r1,nk+1,q1*q2)
        return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)
    return d


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




