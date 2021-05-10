function dmean=DOTNB_mean(r1,p1,r2,p2)

% calculate the variance of the difference of two NB distributions (r1,p1) and (r2,p2)

q1 = 1-p1; q2 = 1-p2;
dmean = vpa(r1*q1/p1^2 + r2*q2/p2^2);

end


