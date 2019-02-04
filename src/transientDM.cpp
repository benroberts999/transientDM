#include "ChronoTimer.h"
#include "ClockNetwork.h"
#include "DataIO.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

enum class WhichOutput { limit, R, limit_and_R, Rthresh };

// void defineIntegerLogGrid(std::vector<int> &grid, int min, int max, int N);
void defineDoubleLogGrid(std::vector<double> &grid, double min, double max,
                         int N);

void dmSearch_tau_int(const ClockNetwork &net, int teff_min, int teff_max,
                      int nteff, TDProfile profile, int nJeffW, double n_sig,
                      double iday, double fday, int min_N_pairs,
                      bool force_PTB_SrYb, bool force_SyrHbNplYb,
                      std::string olabel, WhichOutput whichoutput);

void outputConstraints(const std::vector<double> &jeff_grid, int tau_avg,
                       const std::vector<double> &da_max,
                       const std::vector<double> &Del_da,
                       const std::vector<double> &T_obs, double n_sig,
                       const std::string &olabel);

void outputR_teff(const std::vector<double> &jeff_grid, int tau_avg,
                  const std::vector<double> &R_max,
                  const std::vector<double> &da_max,
                  const std::vector<double> &Del_da,
                  const std::vector<double> &T_obs, const std::string &olabel);

void calculateThreshold(ClockNetwork &net, int teff_min, int teff_max,
                        int nteff, TDProfile profile, int nJeffW, double iday,
                        double fday, int min_N_pairs, bool force_PTB_SrYb,
                        bool force_SyrHbNplYb, std::string olabel);

//******************************************************************************
int main() {
  ChronoTimer timer;

  // Input parameters (to be read from file)
  std::string clock_list_infn;           // File that contains clock list
  double iday, fday;                     // First/last day
  int tau_avg, max_bad_avg;              // average data?
  int teff_min, teff_max, nteff;         // tau_eff grid
  int min_N_pairs;                       // min. allowable # of clock pairs
  bool force_PTB_SrYb, force_SyrHbNplYb; // Force specific clocks?
  TDProfile profile;                     // TD profile (flat/gaussian)
  int n_sig;                             //[1=68%CL, 1.645=90%, 2=95.5%]
  std::string olabel;
  WhichOutput whichoutput; // enum for what output/calculation to do

  // parameters for injecting fake event
  bool inject_fake;
  bool random_inject;
  double da0_inject;
  double tau_inject;
  double day_inject;

  int nJeffW = 1; // window multiplier (not read from file)

  // Read in from input file:
  {
    int iPSY, iSHNY, iP, iwhat, irand;
    auto tp = std::forward_as_tuple(
        clock_list_infn, iday, fday, tau_avg, max_bad_avg, teff_min, teff_max,
        nteff, min_N_pairs, iPSY, iSHNY, iP, n_sig, iwhat, olabel, da0_inject,
        tau_inject, day_inject, irand);
    DataIO::setInputParameters("transientDM.in", tp);
    // Parse some of the input options into ``pc''-friendly form
    force_PTB_SrYb = (iPSY == 1) ? true : false;
    force_SyrHbNplYb = (iSHNY == 1) ? true : false;
    profile = (iP == 1) ? TDProfile::Flat : TDProfile::Gaussian;
    olabel = (olabel == "na") ? "" : "-" + olabel;
    if (iwhat == 0)
      whichoutput = WhichOutput::limit;
    else if (iwhat == 1)
      whichoutput = WhichOutput::R;
    else if (iwhat == 2)
      whichoutput = WhichOutput::limit_and_R;
    else if (iwhat == 3)
      whichoutput = WhichOutput::Rthresh;
    else
      whichoutput = WhichOutput::limit; // default..
    inject_fake = (da0_inject == 0) ? false : true;
    random_inject = (irand == 1) ? true : false;
  }

  // Some checks for valid input data:
  if (teff_min < tau_avg) {
    std::cerr << "teff_min<tau_avg\n";
    return 1;
  }
  if (teff_max == teff_min)
    nteff = 1;
  if (nteff == 1)
    teff_max = teff_min;
  if (min_N_pairs <= 0)
    min_N_pairs = 1;

  // Read-in the list of clock-data filenames
  std::vector<std::string> filenames;
  DataIO::getFileNames(filenames, clock_list_infn);

  // read in the clock data + instantiate ClockNetwork object:
  timer.start();
  ClockNetwork net(filenames, tau_avg, max_bad_avg);
  std::cout << "Time to read files: " << timer.lap_reading_str() << "\n\n";

  bool just_write_data_out = false; // hard-coded..used v. rarey
  if (just_write_data_out) {
    net.writeOutClockData();
    return 0;
  }

  int N_tot_pairs = net.get_NtotPairs();
  if (N_tot_pairs < min_N_pairs) {
    std::cerr << "No clocks? Check data file paths\n";
    return 1;
  }

  // Here: If necisary, inject a fake event (for testing!)
  if (inject_fake) {
    if (random_inject)
      net.replaceWithRandomNoise(FillGaps::no);
    int jeff_inject = (int)(tau_inject / tau_avg);
    long j_earliest = (long)((day_inject)*24 * 60 * 60 / tau_avg);
    long j_latest = (long)((100) * 24 * 60 * 60 / tau_avg);
    std::vector<std::vector<double>> s;
    net.genSignalTemplate(s, nJeffW, jeff_inject, profile);
    int Jw = (int)s[0].size(); // store window
    std::vector<int> indep_pairs;

    long j_in = 0;
    int N_pairs = 0;
    for (long j = j_earliest; j < j_latest; j++) {
      j_in = j;
      net.formIndependentSubnet(indep_pairs, j_in, Jw, 0, force_PTB_SrYb,
                                force_SyrHbNplYb);
      N_pairs = (int)indep_pairs.size();
      if (N_pairs >= min_N_pairs)
        break; // found spot!
    }
    if (N_pairs < min_N_pairs) {
      std::cerr << "Couldn't inject event!\n";
      return 1;
    }
    net.injectFakeEvent(indep_pairs, da0_inject, s, j_in);
    std::cout << "**********************************\n";
    std::cout << "* Injecting a fake event, with:  \n"
              << "*  da0 = " << da0_inject << ", tau_int = " << tau_inject
              << "s,    \n"
              << "*  at day=" << (double)j_in * tau_avg / (24 * 60 * 60)
              << "= epoch: " << j_in << "  \n";
    if (random_inject)
      std::cout << "*  (Injecting into random data). \n";
    else
      std::cout << "*  (Injecting into real data).   \n";
    std::cout << "**********************************\n\n";
  }

  // Output some info to screen:
  std::cout << "Running for days (since MJD:" << CNconsts::MJD_DAY_ZERO
            << "): " << iday << " -> " << fday << " = " << iday * 24 * 60 * 60
            << " -> " << fday * 24 * 60 * 60 << " s.\n";
  std::cout << "With: tau_0 = tau_avg = " << tau_avg << "s.\n";
  std::cout << "There are a total of: " << N_tot_pairs << " clock pairs. "
            << " Only running analysis for N>" << min_N_pairs << "\n";
  std::cout << "Searching between tau_int = " << teff_min << "->" << teff_max
            << "s, in: " << nteff << " steps.\n";
  if (profile == TDProfile::Gaussian)
    std::cout << "Assuming a Gaussian profile.\n";
  else
    std::cout << "Assuming a flat profile.\n";
  std::cout << "\n";

  timer.start();
  switch (whichoutput) {
  case WhichOutput::Rthresh:
    calculateThreshold(net, teff_min, teff_max, nteff, profile, nJeffW, iday,
                       fday, min_N_pairs, force_PTB_SrYb, force_SyrHbNplYb,
                       olabel);
    break;
  default:
    dmSearch_tau_int(net, teff_min, teff_max, nteff, profile, nJeffW, n_sig,
                     iday, fday, min_N_pairs, force_PTB_SrYb, force_SyrHbNplYb,
                     olabel, whichoutput); // yay! So few parameters!
  }
  std::cout << "Time to analyse data: " << timer.lap_reading_str() << "\n";

  std::cout << "\nTotal time: " << timer.reading_str() << "\n";
  return 0;
} // END main()
//******************************************************************************

//******************************************************************************
void defineDoubleLogGrid(std::vector<double> &grid, double min, double max,
                         int N)
/*
Forms a logarithmically-spaced grid of doubles, between [min,max] in N steps.
*/
{
  grid.clear();
  grid.reserve(N);
  grid.push_back(min);
  for (int i = 1; i < N; i++) {
    double x = double(i) / (N - 1);
    grid.emplace_back(min * pow(max / min, x));
  }
}

//******************************************************************************
void dmSearch_tau_int(const ClockNetwork &net, int teff_min, int teff_max,
                      int nteff, TDProfile profile, int nJeffW, double n_sig,
                      double iday, double fday, int min_N_pairs,
                      bool force_PTB_SrYb, bool force_SyrHbNplYb,
                      std::string olabel, WhichOutput whichoutput)
/*
Main routine. Maximises \delta\alpha_0 to find limit, and/or R for search.
*/
{

  bool outputLimit = true;
  bool outputR = false;
  switch (whichoutput) {
  case WhichOutput::limit:
    outputLimit = true;
    outputR = false;
    break;
  case WhichOutput::R:
    outputLimit = false;
    outputR = true;
    break;
  case WhichOutput::limit_and_R:
    outputLimit = true;
    outputR = true;
    break;
  case WhichOutput::Rthresh:
    outputLimit = false;
    outputR = false;
    std::cout << "Rthresh Not implemented yet!\n";
    return;
  }

  int tau_avg = net.get_tau0();
  int max_bad_inJw = 0; // just leave as zero?

  // initial/final epochs to search. (Epochs are time/tau_0.)
  // Convert days -> seconds -> epochs
  long j_init = (long)(iday * 24 * 60 * 60 / tau_avg); // in epochs!
  long j_fin = (long)((fday * 24 * 60 * 60) / tau_avg) + 1;

  // Convert tau_eff (=tau_int) to j_eff (epoch)
  double jeff_min = teff_min / tau_avg;
  double jeff_max = teff_max / tau_avg;
  // Define the tau_eff grid (exponentially spaced grid):
  std::vector<double> jeff_grid;
  // defineIntegerLogGrid(jeff_grid, jeff_min, jeff_max, nteff);
  defineDoubleLogGrid(jeff_grid, jeff_min, jeff_max, nteff);

  // Arrays to store results:
  std::vector<double> da_max(nteff);
  std::vector<double> Del_da(nteff);
  std::vector<double> R_max(nteff);
  std::vector<double> T_obs(nteff);

  int N_tot_pairs = net.get_NtotPairs();

#pragma omp parallel for
  for (size_t it = 0; it < jeff_grid.size(); it++) {
    double j_eff = jeff_grid[it];

    // signal template:
    std::vector<std::vector<double>> s;
    s.reserve(N_tot_pairs);
    net.genSignalTemplate(s, nJeffW, j_eff, profile);
    int Jw = (int)s[0].size(); // store window

    // Loop over j0 (t0): note: jbeg is _beginning_ of Jw window!
    long j_step = 1;     // or tau_eff? tau_eff/2? Makes little diff!
    long num_j_used = 0; // count
    for (long jbeg = j_init; jbeg <= j_fin; jbeg += j_step) {

      // Form the independent sub-network:
      std::vector<int> indep_pairs;
      net.formIndependentSubnet(indep_pairs, jbeg, Jw, max_bad_inJw,
                                force_PTB_SrYb, force_SyrHbNplYb);
      int N_pairs = (int)indep_pairs.size();
      if (N_pairs < min_N_pairs)
        continue;

      ++num_j_used; // count 'good'/used data points (for Tobs)

      auto xHs = net.calculate_dHs_sHs(indep_pairs, s, jbeg);

      // Maximise da and R (sepperately!)
      double da_bf = fabs(xHs.dHs) / xHs.sHs;
      double Delta = 1. / sqrt(xHs.sHs);
      double R = fabs(xHs.dHs) / sqrt(2 * xHs.sHs);
      if (da_bf + Delta > da_max[it] + Del_da[it])
      // if(da_bf>da_max[it])
      {
        da_max[it] = da_bf;
        Del_da[it] = Delta; // 1./sqrt(xHs.sHs);
      }
      if (R > R_max[it])
        R_max[it] = R;

    } // j (t0) loop

    // store observation time (for this tau_eff)
    double Tobs_hr = double(num_j_used * j_step * tau_avg) / (60 * 60.);
    T_obs[it] = Tobs_hr;

  } // tau_eff

  // for tau_eff = 60,120,1020,10020
  // => j_eff = 1,2,17,167
  std::cout << "Summary of results:\n";
  std::cout << "tau_eff  da_max    Da         Rmax   Tobs\n";
  for (size_t i = 0; i < jeff_grid.size() - 1; i++) {
    auto jf = jeff_grid[i];
    auto tf = jf * tau_avg;
    auto tfp = jeff_grid[i + 1] * tau_avg;
    if ((tf <= 60. && tfp > 60.) || (tf <= 100. && tfp > 100.) ||
        (tf <= 1000. && tfp > 1000.) || (tf <= 9900. && tfp > 9900.))
      printf("%5.0fs   %7.1e   %7.1e  %6.2f   %5.2fhr\n", jf * tau_avg,
             da_max[i], Del_da[i], R_max[i], T_obs[i]);
  }
  std::cout << "\n";

  if (outputLimit)
    outputConstraints(jeff_grid, tau_avg, da_max, Del_da, T_obs, n_sig, olabel);

  if (outputR)
    outputR_teff(jeff_grid, tau_avg, R_max, da_max, Del_da, T_obs, olabel);
}

//******************************************************************************
void outputConstraints(const std::vector<double> &jeff_grid, int tau_avg,
                       const std::vector<double> &da_max,
                       const std::vector<double> &Del_da,
                       const std::vector<double> &T_obs, double n_sig,
                       const std::string &olabel) {
  // Define output grids for limits.
  //(hard-coded output grids!)
  int odim_T = 128;
  double oT_min = 0.5; // hours
  double oT_max = 50.; // hours
  std::vector<double> T_grid;
  defineDoubleLogGrid(T_grid, oT_min, oT_max, odim_T);
  int odim_tint = 2048;
  double otint_min = 1.;   // seconds
  double otint_max = 1.e5; // seconds
  std::vector<double> tint_grid;
  defineDoubleLogGrid(tint_grid, otint_min, otint_max, odim_tint);

  // Vector to store output limits
  std::vector<std::vector<double>> dX(odim_T, std::vector<double>(odim_tint));

  // Poisson factor;
  // require: T_obs > T / pois_fac
  // (or, T_max = T_obs * pois_fac)
  double pois_fac = -1. / log(1. - erf(0.7 * n_sig));

  // Sort the limits from the tau_eff grid into the tau_int grid
  // Taking observation time into account for T
  // Note: tau_int << T (otherwise, not transient! take 2x?)
  for (size_t iT = 0; iT < T_grid.size(); iT++) {
    double T = T_grid[iT];
    for (size_t iteff = 0; iteff < jeff_grid.size(); iteff++) {
      double t_eff = jeff_grid[iteff] * tau_avg;
      double Tobs = T_obs[iteff];
      if (T > Tobs * pois_fac)
        continue; // could probs beak, but meh
      double da = da_max[iteff] + n_sig * Del_da[iteff];
      for (size_t itint = 0; itint < tint_grid.size(); itint++) {
        double t_int = tint_grid[itint];
        if (t_int < t_eff)
          continue;
        if (t_int > 0.5 * T * 60. * 60.)
          continue; // tau_int << T [secs/hours]
        double dX_prev = dX[iT][itint];
        if (da < dX_prev || dX_prev == 0)
          dX[iT][itint] = da;
      }
    }
  }

  // Write out the constraints to text file (function of T and tau_int)
  std::string ofname = "dAlpha" + olabel + ".txt";
  std::cout << "Writing out results (Function of T, tau_int) to: " << ofname
            << "\n";
  std::ofstream ofs(ofname);
  ofs << "tau_int\\T ";
  ofs.precision(3);
  for (auto i_T = 0ul; i_T < T_grid.size(); i_T++) {
    ofs << T_grid[i_T] << " ";
  }
  ofs << "\n";
  ofs.precision(3);
  for (auto i_tint = 0u; i_tint < tint_grid.size(); i_tint++) {
    double tau_int = tint_grid[i_tint];
    ofs << tau_int << " ";
    for (auto i_T = 0u; i_T < T_grid.size(); i_T++) {
      ofs << dX[i_T][i_tint] << " ";
    }
    ofs << "\n";
  }
  ofs.close();
}

//******************************************************************************
void outputR_teff(const std::vector<double> &jeff_grid, int tau_avg,
                  const std::vector<double> &R_max,
                  const std::vector<double> &da_max,
                  const std::vector<double> &Del_da,
                  const std::vector<double> &T_obs, const std::string &olabel)
/*
Writes out results (R_max, da_max etc.)
as a function of tau_eff [not tau_int].
Each point will have different T_obs! (also outputted)
*/
{
  // Form:  tau_eff | R da Del_da T_obs
  std::string ofname = "R_teff" + olabel + ".txt";
  std::cout << "Writing out results (Function of tau_eff) to: " << ofname
            << "\n";
  std::ofstream ofs(ofname);
  ofs << "tau_eff R damax Delda Tobs\n";
  ofs.precision(4);

  for (size_t i = 0; i < jeff_grid.size(); i++) {
    ofs << jeff_grid[i] * tau_avg << " " << R_max[i] << " " << da_max[i] << " "
        << Del_da[i] << " " << T_obs[i] << "\n";
  }
  ofs.close();
}

//******************************************************************************
void calculateThreshold(ClockNetwork &net, int teff_min, int teff_max,
                        int nteff, TDProfile profile, int nJeffW, double iday,
                        double fday, int min_N_pairs, bool force_PTB_SrYb,
                        bool force_SyrHbNplYb, std::string olabel)
/*
Calculates R (SNR) threshold.
Simulates data, num_trials times.
Takes the average of the num_avg worst R's
num_trials = num_avg*100
==> roughly, gives 99% C.L
*/
{
  int tau_avg = net.get_tau0();
  int max_bad_inJw = 0; // just leave as zero?

  // initial/final epochs to search. (Epochs are time/tau_0.)
  // Convert days -> seconds -> epochs
  long j_init = (long)(iday * 24 * 60 * 60 / tau_avg); // in epochs!
  long j_fin = (long)((fday * 24 * 60 * 60) / tau_avg) + 1;

  // Convert tau_eff (=tau_int) to j_eff (epoch)
  double jeff_min = teff_min / tau_avg;
  double jeff_max = teff_max / tau_avg;
  // Define the tau_eff grid (exponentially spaced grid):
  std::vector<double> jeff_grid;
  defineDoubleLogGrid(jeff_grid, jeff_min, jeff_max, nteff);

  int N_tot_pairs = net.get_NtotPairs();

  int num_avg = 10;
  int num_trials = num_avg * 100; // 5*num_avgd_trials; //integer multiple!

  std::vector<std::vector<double>> R_trials(num_trials,
                                            std::vector<double>(nteff));

  for (int trial = 0; trial < num_trials; trial++) {

    // Arrays to store results:
    std::vector<double> &R_max = R_trials[trial];

    net.replaceWithRandomNoise(FillGaps::no);

#pragma omp parallel for
    for (size_t it = 0; it < jeff_grid.size(); it++) {
      double j_eff = jeff_grid[it];
      // signal template:
      std::vector<std::vector<double>> s;
      s.reserve(N_tot_pairs);
      net.genSignalTemplate(s, nJeffW, j_eff, profile);
      int Jw = (int)s[0].size(); // store window size
      // Loop over j0 (t0): note: jbeg is _beginning_ of Jw window!
      long j_step = 1; // or tau_eff? tau_eff/2? Makes little diff!
      // long num_j_used = 0; //count
      for (long jbeg = j_init; jbeg <= j_fin; jbeg += j_step) {
        // Form the independent sub-network:
        std::vector<int> indep_pairs;
        net.formIndependentSubnet(indep_pairs, jbeg, Jw, max_bad_inJw,
                                  force_PTB_SrYb, force_SyrHbNplYb);
        if ((int)indep_pairs.size() < min_N_pairs)
          continue;
        // Maximise R
        auto xHs = net.calculate_dHs_sHs(indep_pairs, s, jbeg);
        double R = fabs(xHs.dHs) / sqrt(2 * xHs.sHs);
        if (R > R_max[it])
          R_max[it] = R;
      } // j (t0) loop
    }   // tau_eff

  } // Trials

  std::vector<double> R_max(nteff);
  std::vector<double> R_avg(nteff);
  // //Find the largest (+smallest) num_avg values
  // We want to find the largest n values in M array.
  // (1) Push first n values into new array, A.
  // (2) Sort A.
  //  Loop: from n->M
  //  |--> (3) Get next value.
  //  |    (4) Check if it's larger than one of the current values.
  //  |    (5) If so: insert in into A (in correct location),
  //  |____    and drop last element.
  for (int j = 0; j < nteff; j++) {
    double R_tot_av = 0;
    std::vector<double> max_v; //(num_avg);
    max_v.reserve(num_avg + 1);
    for (int t = 0; t < num_avg; t++) {
      max_v.push_back(R_trials[t][j]); // check??
      R_tot_av += R_trials[t][j];
    }
    std::sort(max_v.begin(), max_v.end(), std::greater<double>());
    for (int t = num_avg; t < num_trials; t++) {
      R_tot_av += R_trials[t][j];
      double newv = R_trials[t][j];
      for (int k = 0; k < num_avg; k++) {
        if (newv > max_v[k]) {
          max_v.insert(max_v.begin() + k, newv);
          max_v.pop_back();
          break;
        }
      }
    }
    // average maximum num_avg values (max_v)
    double sum = std::accumulate(max_v.begin(), max_v.end(), 0.0);
    double mean = sum / num_avg;
    R_max[j] = mean;
    R_tot_av /= num_trials;
    R_avg[j] = R_tot_av;
  }

  std::cout << "Summary of threshold calc:\n";
  std::cout
      << "  nb: R_thresh is calculated as the average of the: \n"
      << num_avg
      << "  worst(largest) R_max's [for each given tau], out of a total of: \n"
      << "  " << num_trials << " trails (roghly, 99 percentile).\n";
  std::cout << "|tau_eff    R_mean    R_thresh|\n";
  for (size_t i = 0; i < jeff_grid.size(); i++) {
    auto jf = jeff_grid[i];
    auto tf = jf * tau_avg;
    auto tfp = jeff_grid[i + 1] * tau_avg;
    if ((tf <= 60. && tfp > 60.) || (tf <= 100. && tfp > 100.) ||
        (tf <= 1000. && tfp > 1000.) || (tf <= 9900. && tfp > 9900.))
      printf("|%5.0fs    %7.4f   %7.4f  |\n", jf * tau_avg, R_avg[i], R_max[i]);
  }
  std::cout << "\n";

  // Form:  tau_eff | R da Del_da T_obs
  std::string ofname = "R_thresh" + olabel + ".txt";
  std::cout << "Writing out results (threshold, fn of tau_eff) to: " << ofname
            << "\n";
  std::ofstream ofs(ofname);
  ofs << "tau_eff R_mean R_thresh\n";
  ofs.precision(4);
  for (size_t i = 0; i < jeff_grid.size(); i++) {
    ofs << jeff_grid[i] * tau_avg << " " << R_avg[i] << " " << R_max[i] << " "
        << "\n";
  }
  ofs.close();
}
