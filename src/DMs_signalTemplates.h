#pragma once

namespace DMsignalTemplate{

double s_Gaussian(double tau0, double tau_int,
  double t0, double tj, double Kab=1.);

double s_topHat(double , double tau_int,
  double t0, double tj, double Kab=1.);

double fastErf(double x);

}//namespace
