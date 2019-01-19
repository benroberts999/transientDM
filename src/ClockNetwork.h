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
    std::vector<long> _total_time;

    std::vector<std::string> _clock_name_A;
    std::vector<std::string> _clock_name_B;

    std::vector<double> _K_A;
    std::vector<double> _K_B;
    std::vector<double> _K_AB;

    // std::vector<double> sigma_tau0;

    //long unsigned int initial_epoch; //?
    //int tau_0;

    // std::vector<bool> pair_OK;

    // std::vector<int> ranked_index_list;



  //private:

    int readInDataFile(const std::string &in_fname);
    void fetchClockInfo(const std::string &fn);

};
