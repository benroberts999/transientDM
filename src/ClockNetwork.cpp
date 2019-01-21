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



// namespace ClockNetwork{
//
// }



//******************************************************************************
int ClockNetwork::get_NtotPairs() const
{
  return (int) _mean.size();
}

//******************************************************************************
ClockNetwork::ClockNetwork(const std::vector<std::string> &filenames,
  int tau_avg, int max_bad)
  : _tau_0(tau_avg)
{

  for(size_t i=0; i<filenames.size(); i++){
    std::cout<<"Reading input data file: "<<i+1<<"/"<<filenames.size()
    // <<" ("<<filenames[i]<<")"
    <<"         \r"<<std::flush;
    readInDataFile(filenames[i],max_bad);
  }
  std::cout<<"\n\n";


  calculateSigma0();
  rankClockPairs();

  std::cout<<"Network summary:\n";
  for(auto i : _ranked_index_list){
    printf("%17s: dK=%5.2f, x0=%8.1e, sig=%8.1e\n",name(i).c_str(),_K_AB[i],
      _mean[i],_sigma0[i]);
  }
  std::cout<<"\n";

  //Reserve space in vectors (this doesn't help too much)
  int expected_num_clocks = (int) filenames.size();
  _delta_omega.reserve(expected_num_clocks);
  _data_ok.reserve(expected_num_clocks);
  _initial_time.reserve(expected_num_clocks);
  _initial_epoch.reserve(expected_num_clocks);
  _mean.reserve(expected_num_clocks);
  _sigma0.reserve(expected_num_clocks);
  _ranked_index_list.reserve(expected_num_clocks);
  _clock_name_A.reserve(expected_num_clocks);
  _clock_name_B.reserve(expected_num_clocks);
  _K_A.reserve(expected_num_clocks);
  _K_B.reserve(expected_num_clocks);
  _K_AB.reserve(expected_num_clocks);

}

//******************************************************************************
std::string ClockNetwork::name(int i)const{
  if(i>=0 && i<(int)_clock_name_A.size())
    return _clock_name_A[i]+"-"+_clock_name_B[i];
  else return "";
}


//******************************************************************************
void ClockNetwork::formIndependentSubnet(std::vector<int> &indep_pairs,
  long beg_epoch, int Jw, int max_bad,
  bool force_PTB_SrYb, bool force_SyrHbNplYb) const
/*
XXX Add ability to force to use PTB-Sr-Yb
*/
{

  std::vector<std::string> clock_list;
  clock_list.reserve(2*_ranked_index_list.size());
  indep_pairs.clear();

  for(auto i : _ranked_index_list){
    //Check bounds:
    long j_beg_i = beg_epoch - _initial_epoch[i];
    if(j_beg_i < 0) continue;//BAD

    long final_epoch_i = _initial_epoch[i] + _delta_omega[i].size();
    if(final_epoch_i < beg_epoch + Jw) continue;
    //XXX Check for off-by-one errors!

    int bad = 0;
    for(long j = 0; j<Jw; j++){
      if(!_data_ok[i][j_beg_i + j]){
        ++bad;
        if(bad > max_bad) break;
      }
    }
    if(bad > max_bad) continue;

    // bool already = (std::find(clock_list.begin(), clock_list.end(),
    // this_clock) != clock_list.end());
    bool already = false;
    for(auto clk : clock_list){
      if (clk == _clock_name_A[i] || clk == _clock_name_B[i]){
        already = true;
        break;
      }
    }
    if(already) continue;

    clock_list.push_back(_clock_name_A[i]);
    clock_list.push_back(_clock_name_B[i]);
    indep_pairs.push_back(i);
  }

  // bool force_PTB_SrYb = true;
  if(force_PTB_SrYb){
    bool found = false;
    for(auto i : indep_pairs){
      if(name(i) == "PTBSr-PTBYb") found = true;
    }
    if(!found) indep_pairs.clear();
  }
  if(force_SyrHbNplYb){
    bool found = false;
    for(auto i : indep_pairs){
      if(name(i) == "SYRTEHg-NPLYb") found = true;
    }
    if(!found) indep_pairs.clear();
  }

}



//******************************************************************************
Result_xHs ClockNetwork::calculate_dHs_sHs(
  const std::vector<int> &indep_pairs,
  const std::vector<std::vector<double> > &s, long beg_epoch) const
/*
beg_epoch is the begining epoch for the window
(epoch is # _points_ since MJD_DAY_ZERO, i.e. time/tau_0)
*/
{
  double dHs = 0.;
  double sHs = 0.;
  auto Jw = s[0].size();
  for(auto i : indep_pairs){
    double Hii = 1./pow(_sigma0[i],2);
    double ss_i=0;
    double ds_i=0;
    long j_beg_i = beg_epoch - _initial_epoch[i];
    for(auto j=0ul; j<Jw; j++){
      double sij = s[i][j];
      double dij = _delta_omega[i][j_beg_i + j];
      ss_i += sij*sij;
      ds_i += dij*sij;
    }
    sHs += Hii*ss_i;
    dHs += Hii*ds_i;
  }
  return Result_xHs(dHs,sHs);
}


//******************************************************************************
void ClockNetwork::genSignalTemplate(std::vector<std::vector<double> > &s,
  int n_window, int j_int, TDProfile profile) const
/*
j_int := tau_int / tau_0
*/
{
  int Jw = n_window * j_int;
  if(Jw%2 == 0) ++Jw; // Jw must be odd

  int tint = j_int * _tau_0;

  s.resize(_K_AB.size(), std::vector<double>(Jw));
  //XXX note: when creating s, _reserve_ enough space for largest JW!!

  double j0 = 0.5*Jw - 1.;
  double t0 = j0*_tau_0;

  //Define a function pointer, to switch between
  double (*dmSignal)(double, double, double, double, double);
  switch(profile){
    case TDProfile::Gaussian : dmSignal = &DMsignalTemplate::s_Gaussian;
      break;
    case TDProfile::Flat : dmSignal = &DMsignalTemplate::s_topHat;
  }


  std::vector<double> sonK;
  sonK.reserve(Jw);
  for(int j=0; j<Jw; j++){
    double tj = _tau_0*j;
    sonK.emplace_back(dmSignal(_tau_0,tint,t0,tj,1));
  }

  for(size_t i=0; i<_K_AB.size(); i++){//loop through clocks
    double Kabi = _K_AB[i];
    for(int j = 0; j<Jw; j++) s[i][j] = Kabi*sonK[j];
  }

}




//******************************************************************************
int ClockNetwork::readInDataFile(const std::string &in_fname, int max_bad)
//XXX Update for other formats!
{

  std::vector<double> times;
  std::vector<double> tmp_dw;
  int ok = DataIO::read_text_XY(in_fname,times,tmp_dw);
  if(ok!=0) return ok;

  int tau_avg = _tau_0;

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

  int mod_avg = (int) (i_time%((long)tau_avg));
  int ibeg = mod_avg==0 ? 0 : tau_avg - mod_avg;

  long new_initial_time = i_time + ibeg;
  if(new_initial_time%tau_avg !=0) std::cerr<<"\nFAIL CN 105\n";

  _initial_time.push_back(new_initial_time);
  _initial_epoch.push_back(new_initial_time/tau_avg);

  //pad-out data with zeros for "missing" points..
  //(and 'mark' bad/missing points)

  //Tranfer data into temporary array, padding missing points w/ zeroes.
  std::vector<double> w_tmp, ok_tmp;
  w_tmp.reserve(tot_time);
  ok_tmp.reserve(tot_time);
  {
    size_t j=0;
    for(long i=0; i<tot_time; i++){
      long tf = (long)round(times[j]);
      long t  = i_time + i;
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

  //std::cout<<"137\n"<<std::flush;

  _delta_omega.push_back({});
  auto &w_ref = _delta_omega.back();
  _data_ok.push_back({});
  auto &ok_ref = _data_ok.back();
  w_ref.reserve(tot_time/tau_avg);
  ok_ref.reserve(tot_time/tau_avg);

  double oa_sum=0;
  long oa_good=0;
  for(long i=ibeg; i<tot_time-tau_avg; i+=tau_avg){
    int bad = 0;
    int good = 0;
    double w_sum = 0;
    bool new_ok = true;
    for(long j=0; j<tau_avg; j++){
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
    }
  }

  //store the mean:
  _mean.emplace_back(oa_sum/double(oa_good));

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

  int i = (int) (fn.find(prefix) + prefix.length());
  int j = (int) fn.find(joiner,i);
  int jl= (int) joiner.length();
  int k = (int) fn.find(suffix,j);

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

  // // Positions: Not used (for now)
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
    _sigma0.push_back(sqrt(var/double(Npts)));
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
  //bool print = true; //Just for testing!

  /*
    * Sorts clocks by
  */
  if(_sigma0.size()==0) std::cerr<<"FAIL CN 268 - no sigma?\n";

  int N_tot_pairs = (int) _sigma0.size();

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
