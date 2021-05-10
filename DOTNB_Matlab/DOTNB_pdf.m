function pvalue=DOTNB_pdf(k,r1,p1,r2,p2)

% calculate the p-value of two different NB distributions (r1,p1) and (r2,p2)
% k is the observed number of difference

q1 = 1-p1; q2 = 1-p2;

if k > 0
    pvalue = p1^r1*p2^r2*((gamma(r1+k)/gamma(r1))/gamma(1+k))*q1^k*hypergeom([r1+k r2],k+1,q1*q2);
else
    k=-k;
    pvalue = p1^r1*p2^r2*((gamma(r2+k)/gamma(r2))/gamma(1+k))*q2^k*hypergeom([r2+k r1],k+1,q1*q2);
end
pvalue=vpa(pvalue);
end


