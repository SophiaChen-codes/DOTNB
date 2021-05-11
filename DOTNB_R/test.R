library(DOTNB)

#the following is to fit a simulated 
r1 = 2
p1 = 0.58
r2 = 3
p2 = 0.3
q1 = 1-p1
q2 = 1-p2

dmean = DOTNB_mean(r1,p1,r2,p2)

dvar = DOTNB_var(r1,p1,r2,p2)

n=30000

diff=rnbinom(n, size=r1, prob=p1)-rnbinom(n, size=r2, prob=p2)

Range=(ceiling(dmean)-30):(ceiling(dmean)+30)
len=length(Range)
prob<- vector(mode='numeric',length=len)

for(i in 1:len){
	prob[i]=DOTNB_pdf(r1,p1,r2,p2,Range[i]);
}

#par(new=TRUE)  #keep the figure, overlay
plot(Range,prob,col="red")
lines(density.default(diff), col = "blue")