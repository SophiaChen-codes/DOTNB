function simulated_array=DOTNB_simulate(num,r1,p1,r2,p2)

% simulate num numbers that generated from two different NB distributions (r1,p1) and (r2,p2)
% num is the number required

q1 = 1-p1; q2 = 1-p2;

simulated_array = zeros(1,num);
for i = 1:num
    simulated_array(i) = nbinrnd(r1,p1)-nbinrnd(r2,p2);
end
simulated_array = vpa (simulated_array);
end


