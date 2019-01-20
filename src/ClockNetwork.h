#pragma once
#include <vector>
#include <string>

const int MJD_DAY_ZERO = 57900;
const int SECS_IN_DAY = 24*60*60;
// class ClockData{
//
// };

enum class ClockAtom{Sr, Hg, YbII};

class ClockNetwork{

  public:

    //ClockNetwork(); //std::vector<std::string> filenames, int tau_avg

    // double getK_AB();

  //private:

    // std::vector<std::vector<double> > raw_data;
    std::vector<std::vector<double> > _delta_omega;
    std::vector<std::vector<bool> > _data_ok;
    std::vector<long> _initial_time;
    std::vector<long> _initial_epoch;
    std::vector<double> _mean;
    std::vector<double> _sigma0;
    // std::vector<long> _total_time;

    std::vector<int> _ranked_index_list;

    int _tau_0 = 1;
    // int tau_avg;

    // std::vector<std::vector<double> > _avgd_delta_omega;
    // std::vector<std::vector<bool> > _avgd_data_ok;

    std::vector<std::string> _clock_name_A;
    std::vector<std::string> _clock_name_B;

    std::vector<double> _K_A;
    std::vector<double> _K_B;
    std::vector<double> _K_AB;

    //int tau0;

    // std::vector<double> sigma_tau0;

    //long unsigned int initial_epoch; //?
    //int tau_0;

    // std::vector<bool> pair_OK;

    // std::vector<int> ranked_index_list;



  //private:

    //int readInDataFile_old(const std::string &in_fname);
    int readInDataFile(const std::string &in_fname,
      int tau_avg=1, int max_bad=0);

    void fetchClockInfo(const std::string &fn);

    void calculateSigma0();
    void rankClockPairs();

    void genSignalTemplate(std::vector<std::vector<double> > &s,
      int n_window, int tint_on_tau0) const;

    std::vector<double> calculate_dHs_sHs(
      const std::vector<int> &indep_pairs,
      const std::vector<std::vector<double> > &s, long beg_epoch) const;

    void formIndependentSubnet(std::vector<int> &indep_pairs,
      long beg_epoch, int Jw, int max_bad) const;

};
