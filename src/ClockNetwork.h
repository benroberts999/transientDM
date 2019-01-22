#pragma once
#include <vector>
#include <string>

namespace CNconsts{
  const int MJD_DAY_ZERO = 57900;
  const int SECS_IN_DAY = 24*60*60;
}

enum class ClockAtom{Sr, Hg, YbII};
enum class TDProfile{Gaussian, Flat};

//******************************************************************************
struct Result_xHs{

  Result_xHs(double dHs, double sHs)
    :dHs(dHs), sHs(sHs){};

  double dHs;
  double sHs;

};

//******************************************************************************
class ClockNetwork{

public:

  ClockNetwork(const std::vector<std::string> &filenames, int tau_avg,
    int max_bad=0);

  //getters:
  int get_NtotPairs() const;
  int get_tau0() const;
  std::string name(int i) const;

  //Used for the analysis:
  void genSignalTemplate(std::vector<std::vector<double> > &s,
    int n_window, int j_int, TDProfile profile) const;

  void formIndependentSubnet(std::vector<int> &indep_pairs,
    long beg_epoch, int Jw, int max_bad,
    bool force_PTB_SrYb=false, bool force_SyrHbNplYb=false) const;

  Result_xHs calculate_dHs_sHs(
    const std::vector<int> &indep_pairs,
    const std::vector<std::vector<double> > &s, long beg_epoch) const;


private: //data

  const int _tau_0; //const means cannot be changed after initialised

  std::vector<std::vector<double> > _delta_omega;
  std::vector<std::vector<bool> > _data_ok;
  std::vector<long> _initial_time; //ever used?
  std::vector<long> _initial_epoch;
  std::vector<double> _mean;
  std::vector<double> _sigma0;

  std::vector<int> _ranked_index_list;

  std::vector<std::string> _clock_name_A;
  std::vector<std::string> _clock_name_B;

  std::vector<double> _K_A;
  std::vector<double> _K_B;
  std::vector<double> _K_AB;


private: //methods

  int readInDataFile(const std::string &in_fname, int max_bad=0);
  void fetchClockInfo(const std::string &fn);
  void calculateSigma0();
  void rankClockPairs();

};
