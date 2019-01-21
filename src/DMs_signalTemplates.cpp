#include "DMs_signalTemplates.h"
#include <cmath>

namespace DMsignalTemplate{

const double sqrtPI = sqrt(M_PI);

//******************************************************************************
double s_Gaussian(double tau0, double tau_int,
  double t0, double tj, double Kab)
/*
\delta\omega/\omega = \varphi(t_j)
\varphi(t_j) = \deltaX * s(t_j)
s(t_j) = 0.5 * K_ab * sqrt(pi) * (tau_int/tau_0)
        * [erf([t0-tj+tau0]/tau_int) - erf([t0-tj]/tau_int)]
The default value of K=1
This is because, often, template is _exactly_ the same for each clock,
only K changes. So only need to call this function for _one_ clock!
*/
{
  if(Kab==0) return 0;
  double amp = 0.5 * Kab * sqrtPI * (tau_int/tau0);
  double arg1 = (t0-tj+tau0)/tau_int;
  double arg2 = (t0-tj)/tau_int;
  return amp*(fastErf(arg1) - fastErf(arg2));
}

//***************************************************************************
double s_topHat(double , double tau_int,
  double t0, double tj, double Kab)
{
  if(Kab==0) return 0;
  double del_t = t0-tj;
  // if (del_t>tau_int) return 0.;
  if(del_t >= -0.5*tau_int && del_t < 0.5*tau_int) return Kab;
  return 0.;
}

//******************************************************************************
double fastErf(double x)
/*
170620.
Up to 10x faster!
Uses a series expansion about 0 for values of |x|<0.35.
Two regions for the series expansion, so it doesn't have to call large powers
in the case that x is very small.
For |x|>0.35, reverts to regular cmath erf function.
Accurate to better than 1e-5
See also:Abramowitz and Stegun (equations 7.1.25–28)
  https://en.wikipedia.org/wiki/Error_function
*/
{
  if(x>3.4) return 1.;
  if(x<-3.4) return -1.;
  if(x<0.005&&x>-0.005) return 1.128379*x;
  if(x<0.2&&x>-0.2) return 1.128379*x-0.3761264*pow(x,3)+0.1128379*pow(x,5);
  if(x<0.35&&x>-0.35)
    return 1.128379*x - 0.3761264*pow(x,3) + 0.1128379*pow(x,5)
    - 0.02686617*pow(x,7);
  double z = fabs(x);
  double w = 0.;
  int sgn = 1;
  if(x<0)sgn = -1;
  if(z<0.7){//order 5 series around 0.5
    double y=z-0.5;
    w=0.5204999+0.8787826*y-0.4393913*pow(y,2)-0.1464638*pow(y,3)
            +0.1830797*pow(y,4)+0.007323188*pow(y,5);
  }else if(z<1.2){//order 5 series around 1
    double y=z-1.;
    w=0.84270+0.415107*y-0.415107*pow(y,2)+0.1383692*pow(y,3)
            +0.0691846*pow(y,4)-0.0691846*pow(y,5);
  }else if(z<2.2){//order 8 series around 1.75
   double y=z-1.75;
    w=0.9866717+0.05277500*y-0.09235624*pow(y,2)+0.09015728*pow(y,3)
            -0.04810221*pow(y,4)+0.006624361*pow(y,5)+0.008963045*pow(y,6)
            -0.006058751*pow(y,7)+0.0007300512*pow(y,8);
  }else if(z<=3.4){//order 7 series around 2.75
    double y=z-2.75;
    w=0.9998994+0.0005862772*y-0.001612262*pow(y,2)+0.002760389*pow(y,3)
            -0.003258114*pow(y,4)+0.002755808*pow(y,5)-0.001657327*pow(y,6)
            +0.0006460410*pow(y,7);
  }
  return sgn*w;
}

}//namespace
