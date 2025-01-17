#include "ClockNetwork.h"
#include "ClockInfo.h"
#include "DMs_signalTemplates.h"
#include "DataIO.h"
#include "RNG_randomNumberGenerators.h"
#include <algorithm> //for sort
#include <string>
#include <vector>

//******************************************************************************
ClockNetwork::ClockNetwork(const std::vector<std::string> &filenames,
                           std::string prefix, std::string suffix,
                           bool skip_bad_bits, int tau_avg, int max_bad)
    : m_tau_0(tau_avg)
// Note: Everywhere inside assumes white noise.
// Therefore, tau_avg must be chosen to be large enough for this to be true!
// There are no checks for this
{

  // Read in and average the raw data.
  for (std::size_t i = 0; i < filenames.size(); i++) {
    std::cout << "Reading input data file: " << i + 1 << "/" << filenames.size()
              << "         \r" << std::flush;
    readInDataFile(filenames[i], max_bad, prefix, suffix, skip_bad_bits);
  }
  std::cout << "\n\n";

  // calculate sigma for each clock (at given tau_avg!) + ranks the clocks
  calculateSigma0();
  rankClockPairs();

  // print a summary of clocks
  std::cout << "Network summary:\n";
  for (int i = 0; i < 61; i++)
    std::cout << "_";
  std::cout << "\n";
  for (auto i : m_ranked_index_list) {
    printf("|%21s: dK=%5.2f, x0=%8.1e, sig=%8.1e |\n", name(i).c_str(),
           m_K_AB[i], m_mean[i], m_sigma0[i]);
  }
  std::cout << "|";
  for (int i = 0; i < 59; i++)
    std::cout << "_";
  std::cout << "|\n";

  // Reserve space in vectors
  // This is only for efficiency, makes no difference to functionality
  // Also, doesn't actually help too much
  auto expected_num_clocks = filenames.size();
  m_delta_omega.reserve(expected_num_clocks);
  m_data_ok.reserve(expected_num_clocks);
  m_initial_time.reserve(expected_num_clocks);
  m_initial_epoch.reserve(expected_num_clocks);
  m_mean.reserve(expected_num_clocks);
  m_sigma0.reserve(expected_num_clocks);
  m_ranked_index_list.reserve(expected_num_clocks);
  m_clock_name_A.reserve(expected_num_clocks);
  m_clock_name_B.reserve(expected_num_clocks);
  m_K_A.reserve(expected_num_clocks);
  m_K_B.reserve(expected_num_clocks);
  m_K_AB.reserve(expected_num_clocks);
}

//******************************************************************************
std::string ClockNetwork::name(int i) const {
  if (i >= 0 && i < (int)m_clock_name_A.size())
    return m_clock_name_A[i] + "-" + m_clock_name_B[i];
  else
    return "";
}

//******************************************************************************
int ClockNetwork::get_NtotPairs() const { return (int)m_mean.size(); }

//******************************************************************************
int ClockNetwork::get_tau0() const { return m_tau_0; }

//******************************************************************************
void ClockNetwork::formIndependentSubnet(std::vector<int> &indep_pairs,
                                         long beg_epoch, int Jw, int max_bad,
                                         bool force_PTB_SrYb,
                                         bool force_SyrHbNplYb) const
// Forms a "sub-network" of only indepedent clock pairs.
// (That is, no single clock is included more than once).
// Stores the _index_'s for each clock-pair in this subnet in indep_pairs
// Note: uses m_ranked_index_list - i.e. which contains a list of _all_
// clock-pair indexes, in order of best to worst! (by sigma^2/kappa)
{

  std::vector<std::string> clock_list;
  clock_list.reserve(2 * m_ranked_index_list.size());
  indep_pairs.clear();
  indep_pairs.reserve(2); // in theory, should be bigger. But i know only ever 2

  for (auto i : m_ranked_index_list) {
    // Check bounds (make sure this clock pair has recorded data during window):
    long j_beg_i = beg_epoch - m_initial_epoch[i];
    if (j_beg_i < 0)
      continue;
    long final_epoch_i = m_initial_epoch[i] + m_delta_omega[i].size();
    if (final_epoch_i < beg_epoch + Jw)
      continue;

    // Count number of 'skipped' points in this window.
    int bad = 0;
    for (long j = 0; j < Jw; j++) {
      if (!m_data_ok[i][j_beg_i + j]) {
        ++bad;
        if (bad > max_bad)
          break;
      }
    }
    if (bad > max_bad)
      continue;

    // Check if one of the two clocks for this pair is already chosen:
    bool already = false;
    for (auto clk : clock_list) {
      if (clk == m_clock_name_A[i] || clk == m_clock_name_B[i]) {
        already = true;
        break;
      }
    }
    if (already)
      continue;

    // If got here, this clock pair is good! Store names to check next clock
    // and store index in output: indep_pairs
    clock_list.push_back(m_clock_name_A[i]);
    clock_list.push_back(m_clock_name_B[i]);
    indep_pairs.push_back(i);
  }

  // Option to ONLY use data when particular clocks are present.
  // This is hard-coded for now....
  // If the required clocks are not present, clear the list. Use no clocks.
  //(there is a faster way to do this, but it's v. messy)
  if (force_PTB_SrYb) {
    bool found = false;
    for (auto i : indep_pairs) {
      if (name(i) == "PTBSr-PTBYb")
        found = true;
    }
    if (!found)
      indep_pairs.clear();
  }
  if (force_SyrHbNplYb) {
    bool found = false;
    for (auto i : indep_pairs) {
      if (name(i) == "SYRTEHg-NPLYb")
        found = true;
    }
    if (!found)
      indep_pairs.clear();
  }
}

//******************************************************************************
void ClockNetwork::writeOutClockData(std::string label) {

  for (auto i : m_ranked_index_list) {
    std::vector<double> x, y;
    auto N = m_delta_omega[i].size();
    x.reserve(N);
    y.reserve(N);
    auto j0 = m_initial_epoch[i];
    for (std::size_t j = 0; j < N; j++) {
      double day = double((j0 + j) * m_tau_0) / (24. * 60 * 60);
      if (!m_data_ok[i][j])
        continue;
      x.push_back(day);
      y.push_back(m_delta_omega[i][j]);
    }
    std::string out_fname = "datplot_" + name(i) + label + ".txt";
    DataIO::write_text_XY(out_fname, x, y);
  }
}

//******************************************************************************
void ClockNetwork::injectFakeEvent(std::vector<int> &indep_pairs, double da0,
                                   std::vector<std::vector<double>> &s,
                                   long beg_epoch) {

  auto Jw = s[0].size();
  for (auto i : indep_pairs) {
    long j_beg_i = beg_epoch - m_initial_epoch[i];
    for (auto j = 0ul; j < Jw; j++) {
      m_delta_omega[i][j_beg_i + j] += da0 * s[i][j];
    }
  }
}

//******************************************************************************
Result_xHs
ClockNetwork::calculate_dHs_sHs(const std::vector<int> &indep_pairs,
                                const std::vector<std::vector<double>> &s,
                                long beg_epoch) const
// beg_epoch is the begining epoch for the window
// (epoch is # _points_ since MJD_DAY_ZERO, i.e. time/tau_0)
// Note: stores the result in a custom struct: "Result_xHs" [defn in CN header]
// XXX Note: assumes diagonal covariance! Must be updated in order to use
// non-diagonal covariance
{
  double dHs = 0.;
  double sHs = 0.;
  auto Jw = s[0].size();
  for (auto i : indep_pairs) {
    double Hii = 1. / pow(m_sigma0[i], 2);
    double ss_i = 0;
    double ds_i = 0;
    long j_beg_i = beg_epoch - m_initial_epoch[i];
    for (auto j = 0ul; j < Jw; j++) {
      double sij = s[i][j];
      double dij = m_delta_omega[i][j_beg_i + j];
      ss_i += sij * sij;
      ds_i += dij * sij;
    }
    sHs += Hii * ss_i;
    dHs += Hii * ds_i;
  }
  return Result_xHs(dHs, sHs); // order here matters: dHs, THEN sHs
}

//******************************************************************************
void ClockNetwork::genSignalTemplate(std::vector<std::vector<double>> &s,
                                     int n_window, double j_int,
                                     TDProfile profile) const
// Creates a signal-template vector, s, for each clock pair:
// Defined:
//   \delta\omega/\omega = \varphi(t_j)
//   \varphi(t_j) = \deltaX * s(t_j)
// Exact form of s depends on profile; given as input (Gaussian or Flat)
// Note: Explicitely assumes coincident case!
// Note: j_int := tau_int / tau_0
{
  // Calculate the window size:
  int Jw = (int)(n_window * j_int);
  // Jw must be odd:
  if (Jw % 2 == 0)
    ++Jw;

  double tint = j_int * m_tau_0;
  double j0 = 0.5 * Jw - 1.;
  double t0 = j0 * m_tau_0;

  // Define a function pointer, to switch between TD profile options
  double (*dmSignal)(double, double, double, double, double);
  dmSignal = &DMsignalTemplate::s_Gaussian;
  if (profile == TDProfile::Flat)
    dmSignal = &DMsignalTemplate::s_topHat;

  // Create vector of s/K_AB.
  // s if different for each clock, but s/K_AB is the same!
  // This is _only_ true for the coincident case
  std::vector<double> sonK;
  sonK.reserve(Jw);
  for (int j = 0; j < Jw; j++) {
    double tj = m_tau_0 * j;
    sonK.emplace_back(dmSignal(m_tau_0, tint, t0, tj, 1.));
  }

  // Write signal-template s
  s.resize(m_K_AB.size(), std::vector<double>(Jw));
  for (std::size_t i = 0; i < m_K_AB.size(); i++) { // loop through clocks
    double Kabi = m_K_AB[i];
    for (int j = 0; j < Jw; j++)
      s[i][j] = Kabi * sonK[j];
  }
}

//******************************************************************************
int ClockNetwork::readInDataFile(const std::string &in_fname, int max_bad,
                                 std::string prefix, std::string suffix,
                                 bool skip_bad_bits)
// Time stored as seconds since MJD.
// Epochs are time / tau_0
{

  std::vector<double> times;
  std::vector<double> tmp_dw;
  int ok =
      DataIO::read_text_XY_conditionalZ(in_fname, times, tmp_dw, skip_bad_bits);
  if (ok != 0)
    return ok;

  int tau_avg = m_tau_0; // already stores in class

  // time stored as seconds since MJD

  // re-scale and shift data points.
  // Input file times in units of days (MJD).
  // I want in units of seconds, since MJD_DAY_ZERO = 57900
  double t_offset = CNconsts::MJD_DAY_ZERO;
  double t_scale = CNconsts::SECS_IN_DAY;
  for (std::size_t i = 0; i < times.size(); i++) {
    times[i] -= t_offset;
    times[i] *= t_scale;
  }

  long i_time = (int)round(times.front());
  long f_time = (int)round(times.back());
  long tot_time = f_time - i_time + 1;

  int mod_avg = (int)(i_time % ((long)tau_avg));
  int ibeg = mod_avg == 0 ? 0 : tau_avg - mod_avg;

  long new_initial_time = i_time + ibeg;
  if (new_initial_time % tau_avg != 0)
    std::cerr << "\nFAIL CN 105\n";

  m_initial_time.push_back(new_initial_time);
  m_initial_epoch.push_back(new_initial_time / tau_avg);

  // pad-out data with zeros for "missing" points..
  //(and 'mark' bad/missing points)
  // Tranfer data into temporary array, padding missing points w/ zeroes.
  std::vector<double> w_tmp;
  std::vector<bool> ok_tmp;
  w_tmp.reserve(tot_time);
  ok_tmp.reserve(tot_time);
  {
    std::size_t j = 0;
    for (long i = 0; i < tot_time; i++) {
      long tf = (long)round(times[j]);
      long t = i_time + i;
      if (tf == t) {
        w_tmp.push_back(tmp_dw[j]);
        ok_tmp.push_back(true);
        j++;
      } else {
        w_tmp.push_back(0.);
        ok_tmp.push_back(false);
      }
    }
  }

  // Copy the temp. arrays into the actual arrays.
  // Note: Technically, faster to all at once (+ avoid large copy)
  // But would be quite complicated
  m_delta_omega.push_back({});
  auto &w_ref = m_delta_omega.back();
  m_data_ok.push_back({});
  auto &ok_ref = m_data_ok.back();
  w_ref.reserve(tot_time / tau_avg);
  ok_ref.reserve(tot_time / tau_avg);
  double oa_sum = 0;
  long oa_good = 0;
  for (long i = ibeg; i < tot_time - tau_avg; i += tau_avg) {
    int bad = 0;
    int good = 0;
    double w_sum = 0;
    bool new_ok = true;
    for (long j = 0; j < tau_avg; j++) {
      if (!ok_tmp[i + j]) {
        ++bad;
        if (bad > max_bad) {
          w_sum = 0.;
          new_ok = false;
          break;
        }
        continue;
      }
      w_sum += w_tmp[i + j];
      ++good;
    }
    double new_w = good > 0 ? w_sum / good : 0;
    w_ref.push_back(new_w);
    ok_ref.push_back(new_ok);
    if (new_ok) {
      oa_sum += new_w;
      ++oa_good;
    }
  }

  // store the mean:
  m_mean.emplace_back(oa_sum / double(oa_good));
  // Get clock info (K's, lab positions etc.)
  fetchClockInfo(in_fname, prefix, suffix);

  return 0;
}

//******************************************************************************
void ClockNetwork::fetchClockInfo(const std::string &fn, std::string prefix,
                                  std::string suffix) {
  std::string joiner = "-";

  int i = (int)(fn.find(prefix) + prefix.length());
  int j = (int)fn.find(joiner, i);
  int jl = (int)joiner.length();
  int k = (int)fn.find(suffix, j);

  std::string clka = fn.substr(i, j - i);
  std::string clkb = fn.substr(j + jl, k - j - jl);

  m_clock_name_A.push_back(clka);
  m_clock_name_B.push_back(clkb);

  int ia = ClockInfo::GET_CLOCK_INDEX(clka);
  int ib = ClockInfo::GET_CLOCK_INDEX(clkb);

  double Ka = ClockInfo::KALPHA[ia];
  double Kb = ClockInfo::KALPHA[ib];

  m_K_A.push_back(Ka);
  m_K_B.push_back(Kb);
  m_K_AB.push_back(Ka - Kb);

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
void ClockNetwork::calculateSigma0() {
  if (m_mean.size() == 0)
    std::cerr << "FAIL CN 237 - no mean?\n";

  m_sigma0.reserve(m_delta_omega.size());
  for (std::size_t i = 0; i < m_delta_omega.size(); i++) {
    double x0 = m_mean[i];
    double var = 0.;
    long Npts = 0;
    for (std::size_t j = 0; j < m_delta_omega[i].size(); j++) {
      if (!m_data_ok[i][j])
        continue;
      var += pow(m_delta_omega[i][j] - x0, 2);
      ++Npts;
    }
    m_sigma0.push_back(sqrt(var / double(Npts)));
  }
}

//******************************************************************************
bool sortcol(const std::vector<double> &v1, const std::vector<double> &v2) {
  // Lambda to sort cols of one vector by another
  return v1[0] > v2[0];
}
//******************************************************************************
void ClockNetwork::rankClockPairs()
// "Sorts" the clock pairs by "Goodness"
// Goodness is defined as K_AB/sigma^2 - effective sensitivity of each clock.
// Note: doesn't sort the actual data, just stores the indexes of each clock
// in the correct 'sorted' order, in the member variable list:
// m_ranked_index_list
{
  if (m_sigma0.size() == 0)
    std::cerr << "FAIL CN 268 - no sigma?\n";

  int N_tot_pairs = (int)m_sigma0.size();

  std::vector<std::vector<double>> m;
  m.reserve(N_tot_pairs);
  for (int i = 0; i < N_tot_pairs; i++) {
    double m_tmp = pow(m_sigma0[i], 1) / fabs(m_K_AB[i]);
    m.emplace_back(std::initializer_list<double>{m_tmp, (double)i + 0.1});
    //+0.1 to prevent rounding error when going from double -> int
  }

  std::sort(m.rbegin(), m.rend(), sortcol);

  for (int i = 0; i < N_tot_pairs; i++)
    m_ranked_index_list.push_back(int(m[i][1]));
}

//******************************************************************************
void ClockNetwork::replaceWithRandomNoise(FillGaps fill_gapsQ)
// Replaces clock data with Random white (Gaussian) noise, with the same
// mean and standard deviation as the existing clock data.
// fill_gapsQ = FillGaps::yes, it will fill-in all the gaps.
// Otherwise, data gaps will be preserved (defult is no)
{

  bool fill_gaps = fill_gapsQ == FillGaps::yes ? true : false;

  if (fill_gaps) {
#pragma omp parallel for
    for (std::size_t i = 0; i < m_delta_omega.size(); i++) {
      double x0 = m_mean[i];
      double sig = m_sigma0[i];
      for (std::size_t j = 0; j < m_delta_omega[i].size(); j++) {
        m_data_ok[i][j] = true;
        m_delta_omega[i][j] = RNG::randGausVal(sig, x0);
      }
    }
  } else {
#pragma omp parallel for
    for (std::size_t i = 0; i < m_delta_omega.size(); i++) {
      double x0 = m_mean[i];
      double sig = m_sigma0[i];
      for (std::size_t j = 0; j < m_delta_omega[i].size(); j++) {
        if (!m_data_ok[i][j])
          continue;
        m_delta_omega[i][j] = RNG::randGausVal(sig, x0);
      }
    }
  }
}
