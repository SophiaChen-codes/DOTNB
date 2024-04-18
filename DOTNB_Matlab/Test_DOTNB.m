%This is to plot simulated numbers and theoretical DOTNB distribution
% parameters used to test
r1 = 2; p1 = 0.58; r2 = 3; p2 = 0.3;
q1 = 1-p1; q2 = 1-p2;

%Calculate the theoretical mean and variance
T_mean=DOTNB_mean(r1,p1,r2,p2)
T_var=DOTNB_var(r1,p1,r2,p2)

%Generate num numbers following DOTNB distribution.
num = 30000
range = double(ceil(T_mean)-30:ceil(T_mean)+30); % plot up and down 30 points arround mean
data = zeros(1,num);
data=double(DOTNB_simulate(num,r1,p1,r2,p2));
[x,y] = hist(data,range);
plot(y,x/sum(x)); hold on

len = length(range); z = zeros(1,len);
% this is used to test a real negative number
test_real=-22.2 
ss= DOTNB_pdf(test_real,r1,p1,r2,p2)
% end of test

%Plot the theoretical pdf
for k = 1:len
    z(k) = DOTNB_pdf(range(k),r1,p1,r2,p2);    
end
 plot(range,z,'r'); hold on   
   
 