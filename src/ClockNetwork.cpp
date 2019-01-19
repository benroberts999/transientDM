#include "ClockNetwork.h"
#include "ClockInfo.h"
#include "DataIO.h"
#include <vector>
#include <string>


//******************************************************************************
int ClockNetwork::readInDataFile(const std::string &in_fname)
//XXX Update for other formats!
{

  std::vector<double> times;
  std::vector<double> tmp_dw;
  std::cout<<"15\n"<<std::flush;
  int ok = DataIO::read_text_XY(in_fname,times,tmp_dw);
  if(ok!=0) return ok;
  std::cout<<"18\n"<<std::flush;

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

  _initial_time.push_back(i_time);
  _total_time.push_back(tot_time);

  //pad-out data with zeros for "missing" points..
  //(and 'mark' bad/missing points)
  _delta_omega.push_back({});
  auto &w_ref = _delta_omega.back();
  _data_ok.push_back({});
  auto &ok_ref = _data_ok.back();
  w_ref.reserve(tot_time);
  ok_ref.reserve(tot_time);
  {
    size_t j=0;
    for(int i=0; i<tot_time; i++){
      int tf = (int)round(times[j]);
      int t  = i_time + i;
      if(tf==t){
        w_ref.push_back(tmp_dw[j]);
        ok_ref.push_back(true);
        j++;
      }else{
        w_ref.push_back(0.);
        ok_ref.push_back(false);
      }
    }
  }
  std::cout<<"64\n"<<std::flush;


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
