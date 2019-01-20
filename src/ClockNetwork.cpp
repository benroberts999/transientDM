#include "ClockNetwork.h"
#include "ClockInfo.h"
#include "DataIO.h"
#include "DMs_signalTemplates.h"
#include <vector>
#include <string>
#include <algorithm>


// //******************************************************************************
// int ClockNetwork::readInDataFile_old(const std::string &in_fname)
// //XXX Update for other formats!
// {
//
//   std::vector<double> times;
//   std::vector<double> tmp_dw;
//   int ok = DataIO::read_text_XY(in_fname,times,tmp_dw);
//   if(ok!=0) return ok;
//
//   //initial_time = (times.front() - MJD_DAY_ZERO) * SECS_IN_DAY;
//
//   //time stored as seconds since MJD
//
//   //re-scale and shift data points.
//   //Input file times in units of days (MJD).
//   //I want in units of seconds, since MJD_DAY_ZERO = 57900
//   double t_offset = MJD_DAY_ZERO;
//   double t_scale = SECS_IN_DAY;
//   for(size_t i=0; i<times.size(); i++){
//     times[i] -= t_offset;
//     times[i] *= t_scale;
//   }
//
//   long i_time = (int)round(times.front());
//   long f_time = (int)round(times.back());
//   long tot_time = f_time - i_time + 1;
//
//   _initial_time.push_back(i_time);
//   _total_time.push_back(tot_time);
//
//   //pad-out data with zeros for "missing" points..
//   //(and 'mark' bad/missing points)
//   _delta_omega.push_back({});
//   auto &w_ref = _delta_omega.back();
//   _data_ok.push_back({});
//   auto &ok_ref = _data_ok.back();
//   w_ref.reserve(tot_time);
//   ok_ref.reserve(tot_time);
//   {
//     size_t j=0;
//     for(int i=0; i<tot_time; i++){
//       int tf = (int)round(times[j]);
//       int t  = i_time + i;
//       if(tf==t){
//         w_ref.push_back(tmp_dw[j]);
//         ok_ref.push_back(true);
//         j++;
//       }else{
//         w_ref.push_back(0.);
//         ok_ref.push_back(false);
//       }
//     }
//   }
//
//   fetchClockInfo(in_fname);
//
//   return 0;
// }


//******************************************************************************
void ClockNetwork::genSignalTemplate(std::vector<std::vector<double> > &s,
  int n_window, int tint_on_tau0) const
{
  int Jw = n_window * tint_on_tau0;
  if(Jw%2 == 0) ++Jw; // Jw must be odd

  int tint = tint_on_tau0 * _tau_0;
  //XXX also have integer tint_on_tau0 ?

  s.resize(_K_AB.size(), std::vector<double>(Jw));

  // int t0 = Tw / 2 / _tau_0;

  double j0 = 0.5*Jw - 1.;
  double t0 = j0*_tau_0;

  std::cout<<"t0="<<t0<<", Jw="<<Jw<<", Tw="<<Jw*_tau_0<<"\n";

  for(size_t i=0; i<_K_AB.size(); i++){//loop through clocks
    double Kab = _K_AB[i];
    for(int j = 0; j < Jw; j++){
      double tj = _tau_0*j;
      //s[i][j] = DMsignalTemplate::s_Gaussian(Kab,_tau_0,tint,t0,tj);
      s[i][j] = DMsignalTemplate::s_topHat(Kab,_tau_0,tint,t0,tj);
    }
  }
}




//******************************************************************************
int ClockNetwork::readInDataFile(const std::string &in_fname,
  int tau_avg, int max_bad)
//XXX Update for other formats!
{

  std::vector<double> times;
  std::vector<double> tmp_dw;
  int ok = DataIO::read_text_XY(in_fname,times,tmp_dw);
  if(ok!=0) return ok;

  _tau_0 = tau_avg; //XXX put this outside? XXX

  //initial_time = (times.front() - MJD_DAY_ZERO) * SECS_IN_DAY;

  //time stored as seconds since MJD

  //re-scale and shift data points.
  //Input file times in units of days (MJD).
  //I want in units of seconds, since MJD_DAY_ZERO = 57900
  double t_offset = MJD_DAY_ZERO;
  double t_scale = SECS_IN_DAY;
  for(size_t i=0; i<times.size(); i++){
    times[i] -= t_offset;
    times[i] *= t_scale;
  }

  long i_time = (int)round(times.front());
  long f_time = (int)round(times.back());
  long tot_time = f_time - i_time + 1;

  //_initial_time.push_back(i_time);
  //_total_time.push_back(tot_time);

  int mod_avg = i_time%tau_avg;
  int ibeg = mod_avg==0 ? 0 : tau_avg - mod_avg;

  long new_initial_time = i_time + ibeg;
  if(new_initial_time%tau_avg !=0) std::cerr<<"\nFAIL CN 105\n";

  _initial_time.push_back(new_initial_time);

  //pad-out data with zeros for "missing" points..
  //(and 'mark' bad/missing points)

  //Tranfer data into temporary array, padding missing points w/ zeroes.
  std::vector<double> w_tmp, ok_tmp;
  w_tmp.reserve(tot_time);
  ok_tmp.reserve(tot_time);
  {
    size_t j=0;
    for(int i=0; i<tot_time; i++){
      int tf = (int)round(times[j]);
      int t  = i_time + i;
      if(tf==t){
        w_tmp.push_back(tmp_dw[j]);
        ok_tmp.push_back(true);
        j++;
      }else{
        w_tmp.push_back(0.);
        ok_tmp.push_back(false);
      }
    }
  }

  std::cout<<"137\n"<<std::flush;

  _delta_omega.push_back({});
  auto &w_ref = _delta_omega.back();
  _data_ok.push_back({});
  auto &ok_ref = _data_ok.back();
  w_ref.reserve(tot_time/tau_avg);
  ok_ref.reserve(tot_time/tau_avg);

  //int max_bad = 0;
  double oa_sum=0;
  long oa_good=0;
  long oa_bad=0;
  for(int i=ibeg; i<tot_time-tau_avg; i+=tau_avg){
    int bad = 0;
    int good = 0;
    double w_sum = 0;
    bool new_ok = true;
    for(int j=0; j<tau_avg; j++){
      if(!ok_tmp[i+j]){
        ++bad;
        if(bad > max_bad){
          w_sum = 0.;
          new_ok = false;
          break;
        }
        continue;
      }
      w_sum += w_tmp[i+j];
      ++good;
    }
    double new_w = good>0? w_sum/good : 0;
    w_ref.push_back(new_w);
    ok_ref.push_back(new_ok);
    if(new_ok){
      oa_sum += new_w;
      ++oa_good;
    }else{
      ++oa_bad;
    }
  }

  double mean = oa_sum/oa_good;
  //std::cout<<mean<<"\n";

  _mean.push_back(mean);

  // std::cout<<"\nhere:\n";
  // std::cout<<oa_good<<" "<<oa_bad<<" = "<<oa_good+oa_bad<<" "<<w_ref.size()
  //   <<" "<<ok_ref.size()<<" |"<<w_ref.size()*tau_avg<<" =? "<<w_tmp.size()<<"\n";

  fetchClockInfo(in_fname);

  return 0;
}

//******************************************************************************
void ClockNetwork::fetchClockInfo(const std::string &fn)
{
  //XXX Update for other formats!
  std::string prefix="cppdiff_";
  std::string joiner="-";
  std::string suffix=".dat";

  int i = fn.find(prefix) + prefix.length();
  int j = fn.find(joiner,i);
  int jl= joiner.length();
  int k = fn.find(suffix,j);

  std::string clka = fn.substr(i,j-i);
  std::string clkb = fn.substr(j+jl,k-j-jl);

  _clock_name_A.push_back(clka);
  _clock_name_B.push_back(clkb);

  int ia = ClockInfo::GET_CLOCK_INDEX(clka);
  int ib = ClockInfo::GET_CLOCK_INDEX(clkb);

  double Ka = ClockInfo::KALPHA[ia];
  double Kb = ClockInfo::KALPHA[ib];

  _K_A.push_back(Ka);
  _K_B.push_back(Kb);
  _K_AB.push_back(Ka-Kb);

  // rA.resize(3);
  // rB.resize(3);
  // double tL = 0;
  // for(int ix=0; ix<3; ix++){
  //   tL += pow(CLK::LABPOS[ia][ix]-CLK::LABPOS[ib][ix],2);
  //   //Absolute position: (relative to SYRTE):
  //   rA[ix] = CLK::LABPOS[ia][ix] - CLK::SYRTE_POS[ix];
  //   rB[ix] = CLK::LABPOS[ib][ix] - CLK::SYRTE_POS[ix];
  // }
  // L = sqrt(tL);

  // tau0 = CLK::tau0;
  // tau  = tau0;

}

//******************************************************************************
void ClockNetwork::calculateSigma0()
{
  if(_mean.size()==0) std::cerr<<"FAIL CN 237 - no mean?\n";

  _sigma0.reserve(_delta_omega.size());
  for(size_t i=0; i<_delta_omega.size(); i++){
    double x0 = _mean[i];
    double var = 0.;
    long Npts = 0;
    for(size_t j=0; j<_delta_omega[i].size(); j++){
      if(!_data_ok[i][j]) continue;
      var += pow(_delta_omega[i][j]-x0,2);
      ++ Npts;
    }
    _sigma0.push_back(sqrt(var/Npts));
  }
}



//******************************************************************************
bool sortcol( const std::vector<double>& v1,
               const std::vector<double>& v2 ) {
    return v1[0] > v2[0];
}
//******************************************************************************
void ClockNetwork::rankClockPairs()
{
  bool print = true; //Just for testing!

  /*
    * Sorts clocks by
  */
  if(_sigma0.size()==0) std::cerr<<"FAIL CN 268 - no sigma?\n";

  int N_tot_pairs = _sigma0.size();

  std::vector< std::vector<double> > m;
  for(int i=0; i<N_tot_pairs; i++){
    double m_tmp = _sigma0[i]/fabs(_K_AB[i]);
    m.push_back({m_tmp,(double)i+0.1});
    //+0.1 to prevent rounding error when going from double -> int
  }

  // if(print){
  //   std::cout<<"pre-sorting:\n";
  //   for(int i=0; i<N_tot_pairs; i++)
  //     std::cout<<m[i][1]<<" "<<m[i][0]<<" "<<_K_AB[i]
  //     <<" "<<_sigma0[i]<<"\n";
  //   std::cout<<"\n";
  // }

  std::sort(m.rbegin(), m.rend(), sortcol);
  for(int i=0; i<N_tot_pairs; i++){
    _ranked_index_list.push_back(int(m[i][1]));
  }

  // if(print){
  //   std::cout<<"post-sorting:\n";
  //   for(int i=0; i<N_tot_pairs; i++)
  //   std::cout<<m[i][1]<<" "<<m[i][0]<<" "<<_K_AB[i]
  //     <<" "<<_sigma0[i]<<"\n";
  //   std::cout<<"\n";
  // }

}
