#include "stdlib.h"
#include "stdio.h"
#include "math.h"



double vs_erfc(double x) 

{
int i;
double sum,unit,x2,s(1.0),rval,l;
double sqpi=sqrt(3.1415926);
x2=x*x;

if(x<-1.99) return 1.99;

if(fabs(x)<2.0) {
  unit=x;
  sum=x;
  double l=2.0/sqpi;
  for(i=1;i<30;i++) {
    s=(-1.0)*s;
    unit=unit*x2/(float)i;
    sum=sum+s*unit/(float)(2*i+1);
    }
  return 1.0-l*sum;
  }
else {
  double l=exp(-x2)/(sqpi*x);
  unit=1.0;
  sum=1.0;
  x2*=2.0;
  for(i=1;i<30;i++) {
    s=(-1.0)*s;
    unit=(float)(2*i-1)*unit/x2;
    sum=sum+s*unit;
    }
  return l*sum;
  }

}
