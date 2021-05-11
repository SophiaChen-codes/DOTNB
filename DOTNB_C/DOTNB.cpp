#include <iostream>
#include <cmath>
using namespace std;

// Note convergence restrictions: abs(x) < 1 and c not a negative integer or zero

double hypergeometric( double a, double b, double c, double x )
{
   const double TOLERANCE = 1.0e-10;
   double term = a * b * x / c;
   double value = 1.0 + term;
   int n = 1;

   while ( abs( term ) > TOLERANCE )
   {
      a++, b++, c++, n++;
      term *= a * b * x / c / n;
      value += term;
   }
   return value;
}

double DOTNB_pdf(int r1, double p1, int r2, double p2, int dn )
{ // calculate the p-value of two different NB distributions (r1,p1) and (r2,p2)
  // dn is the observed number of difference
   double q1 = 1-p1, q2 = 1-p2;
   double prob;
   if(dn>0){
     prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r1+dn)/tgamma(r1))/tgamma(1+dn))*pow(q1,dn)*hypergeometric(r1+dn,r2,dn+1,q1*q2);
   }else{
      dn=-dn;
      prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r2+dn)/tgamma(r2))/tgamma(1+dn))*pow(q2,dn)*hypergeometric(r2+dn,r1,dn+1,q1*q2);
   }
   return prob;
}

//======================================================================

int main()
{
   int r1=2, r2=3,k=10;
   double p1=0.58,p2=0.3;
   cout<<"The probability of the number "<<k<< " of difference occurring is "<< DOTNB_pdf(r1,p1,r2,p2,k)<<"\n";
}
