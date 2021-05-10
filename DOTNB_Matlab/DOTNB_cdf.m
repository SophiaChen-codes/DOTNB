function cvalue=DOTNB_cdf(k,r1,p1,r2,p2)

% calculate the cumulative probability of k that is the difference of two NB distributions (r1,p1) and (r2,p2)
% k is the observed number of difference
%sum pvalue for <= k, 
%use a small stop_pvalue as threshold for stop condition of interation

stop_pvalue=1.0e-20;
cvalue=vpa(0);

if k < DOTNB_mean(r1,p1,r2,p2)
    if DOTNB_pdf(k,r1,p1,r2,p2) < stop_pvalue
        cvalue=cvalue+DOTNB_pdf(k,r1,p1,r2,p2);
    else
        while(DOTNB_pdf(k,r1,p1,r2,p2) > stop_pvalue)
            cvalue=cvalue+DOTNB_pdf(k,r1,p1,r2,p2);
            k=k-1;
        end
    end
else
    if DOTNB_pdf(k,r1,p1,r2,p2) < stop_pvalue
        cvalue=1;
    else
        while(DOTNB_pdf(k,r1,p1,r2,p2) > stop_pvalue)
            cvalue=cvalue+DOTNB_pdf(k,r1,p1,r2,p2);
            k=k-1;
        end
    end
end
   
end


