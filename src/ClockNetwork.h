#pragma once
#include <string>
#include <vector>

namespace CNconsts {
static const int MJD_DAY_ZERO = 57900;
static const int SECS_IN_DAY = 24 * 60 * 60;
} // namespace CNconsts

enum class ClockAtom { Sr, Hg, YbII };
enum class TDProfile { Gaussian, Flat };
enum class FillGaps { yes, no };

/*
This uses two assumptions:
  1) All clocks affected by TD simultaniously (good for tau >~ 15s)
  2) All noise in white Gaussian frequency noise

To correct (1), need to update:
  * Need to store positions of all clocks (this is done, but not used)
  * genSignalTemplate() routine needs to be updated to take into account
  the position of each clock, and the incident speed/direction of TD
  (see my PRD 2016 paper for formulas)

To correct (2), need to update:
 * Need to calculate (or read in) covariance matrix + inverse
 * Need to update calculate_dHs_sHs() function to use matrix

*/

//******************************************************************************
struct Result_xHs {
  // simple struct to store dHs and sHs values
  Result_xHs(double dHs, double sHs) : dHs(dHs), sHs(sHs){};

  double dHs;
  double sHs;
};

//******************************************************************************
class ClockNetwork {

public:
  ClockNetwork(const std::vector<std::string> &filenames, std::string prefix,
               std::string suffix, bool skip_bad_bits, int tau_avg,
               int max_bad = 0);

  // getters:
  int get_NtotPairs() const;
  int get_tau0() const;
  std::string name(int i) const;

  void injectFakeEvent(std::vector<int> &indep_pairs, double da0,
                       std::vector<std::vector<double>> &s, long beg_epoch);

  void replaceWithRandomNoise(FillGaps fill_gapsQ = FillGaps::no);

  // Used for the analysis:
  void genSignalTemplate(std::vector<std::vector<double>> &s, int n_window,
                         double j_int, TDProfile profile) const;

  void formIndependentSubnet(std::vector<int> &indep_pairs, long beg_epoch,
                             int Jw, int max_bad, bool force_PTB_SrYb = false,
                             bool force_SyrHbNplYb = false) const;

  Result_xHs calculate_dHs_sHs(const std::vector<int> &indep_pairs,
                               const std::vector<std::vector<double>> &s,
                               long beg_epoch) const;

  void writeOutClockData(std::string label = "");

private: // data
  const int m_tau_0;

  std::vector<std::vector<double>> m_delta_omega;
  std::vector<std::vector<bool>> m_data_ok;
  std::vector<long> m_initial_time; // ever used?
  std::vector<long> m_initial_epoch;
  std::vector<double> m_mean;
  std::vector<double> m_sigma0;

  std::vector<int> m_ranked_index_list;

  std::vector<std::string> m_clock_name_A;
  std::vector<std::string> m_clock_name_B;

  std::vector<double> m_K_A;
  std::vector<double> m_K_B;
  std::vector<double> m_K_AB;

private: // methods
  int readInDataFile(const std::string &in_fname, int max_bad,
                     std::string prefix, std::string suffix,
                     bool skip_bad_bits);
  void fetchClockInfo(const std::string &fn, std::string prefix,
                      std::string suffix);
  void calculateSigma0();
  void rankClockPairs();
};
