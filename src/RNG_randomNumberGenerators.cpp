#include "RNG_randomNumberGenerators.h"
#include <cmath>
#include <random>
#include <thread>

namespace RNG{

// thread_local std::mt19937 generator(std::random_device{}());

//******************************************************************************
double randDouble(double a, double b)
/*
Returns a uniformly distributed real number between a and b.
Uses 'std::mt19937', a c++11 Mersenne Twister pseudo-random generator.
Also uses <thread>, and is safe to use with OpenMP parallelisation
It does not need to be seeded, that is taken care of automatically.
#include <random>
#include <thread>
*/
{
  thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> distribution(a, b);
  return distribution(generator);
}

//******************************************************************************
double randGausVal(double sig, double x0)
/*
170713.
Uses inverse transform sampling to create a Gaussian random number.
Returns a random double that is drawn from a Gaussian PDF, with mean x0,
and standard deviation 'sig'.
Makes use of my 'randDouble' function.
It also makes use of my "fastInvErf" function.
MUCH faster than previous GSL version (many orders of magnitude)!!
Also, using 'fastInvErf' is ~3x faster than inverseErf.
Output was tested using 10,000 points. Histogram matched exactly! :)
Works with sig=0 (just returns x0)!
=== Change Log ===
*/
{
  double u = randDouble(0,1);  //uniform u
  return x0 + 1.41421*sig*fastInvErf(2*u-1);
}

//******************************************************************************
double fastInvErf(double x)
/*
170622.
Function uses an approximate method to calculate the inverse error function.
There is no cmath inverseErf function.
x must be in interval (-1,1)
It will return '+/-100' if -1 or 1 is given. ok?
It uses a series expansion for small |x|<0.25.
And for 0.25<|x|<0.95, uses series around non-0 points.
Several regions for the series expansion, so it doesn't have to call large
powers in the case that x is very small etc.
Then, uses an algorith I found from:
https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
--I couldn't find the original source
Could be made faster by using an approximation for the log!
This method is ~8x faster than other for |x|<0.95
Accurate to 1e-4 (almost 1e-5)
*/
{
  if(x<0.01&&x>-0.01)return 0.886227*x;
  if(x<0.25&&x>-0.25)return 0.886227*x-0.2320137*pow(x,3)+0.1275562*pow(x,5);
  if(x>=1)return 100.;   //?? allows for error due to floating point errors
  if(x<=-1)return -100.; //??
  double z=fabs(x);
  int sgn=1;
  if(x<0)sgn=-1;
  if(z<=0.95){
    double w;
    if(z<0.55){//order 7 series around 0.4
      double y=z-0.4;
      w=0.3708072+1.016856*y+0.383413*pow(y,2)+0.543234*pow(y,3)
        +0.571544*pow(y,4)+0.777862*pow(y,5)+1.028110*pow(y,6)
        +1.45202*pow(y,7);
    }else if(z<0.85){//order 8 series around 0.7
      double y=z-0.7;
      w=0.7328691+1.516363*y+1.68513*pow(y,2)+3.65912*pow(y,3)+8.6827*pow(y,4)
        +22.4762*pow(y,5)+60.945*pow(y,6)+170.820*pow(y,7)+490.30*pow(y,8);
    }else{//order 8 series around 0.9
      double y=z-0.9;
      w=1.163087+3.42804*y+13.6680*pow(y,2)+86.089*pow(y,3)+621.95*pow(y,4)
        +4846.6*pow(y,5)+39583.*pow(y,6)+3.3382e5*pow(y,7)+2.8817e6*pow(y,8);
    }
    return sgn*w;
  }
  //if z>0.95, just use "normal" approximation.
  double lnx=log(1.-x*x);
  double tt1=4.33+0.5*lnx;
  double tt2=6.803*lnx;
  return sgn*sqrt(sqrt(tt1*tt1-tt2)-tt1);
}

}//namespace
