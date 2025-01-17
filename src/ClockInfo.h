#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace ClockInfo {

const int tau0 = 1; // Intrinsic sampling time from files!
// I assume they're all the same!

const std::vector<std::string> CLOCKS = {
    "NPLSr",    "NPLYb",    "PTBSr",      "PTBYb",     "SYRTEHg",
    "SYRTESr2", "SYRTESrB", "SYRTEFO2Cs", "SYRTEFO2Rb"};

// Absolute frequencies, in Hz [don't actualy use these]
const double NU0_SR = 429228004229873.0;  // Sr
const double NU0_YB = 642121496772645.0;  // Yb
const double NU0_HG = 1128575290808154.4; // Hg
const double NU0_CS = 0.0;                // Cs
const double NU0_RB = 0.0;                // Rb
//
const std::vector<double> NU0 = {NU0_SR, NU0_YB, NU0_SR, NU0_YB, NU0_HG,
                                 NU0_SR, NU0_SR, NU0_CS, NU0_RB};

// Sensitivity coeficients, optical transition. Including Ry.
const double KALPHA_SR = 2. + 0.06;
const double KALPHA_YB = 2. - 6.0;
const double KALPHA_HG = 2. + 0.81;
const double KALPHA_CS = 4. + 0.83;
const double KALPHA_RB = 4. + 0.34;
//
const std::vector<double> KALPHA = {KALPHA_SR, KALPHA_YB, KALPHA_SR,
                                    KALPHA_YB, KALPHA_HG, KALPHA_SR,
                                    KALPHA_SR, KALPHA_CS, KALPHA_RB};

// Lab (clock) positions, in km.
// NOTE: These are not that accurate (city level) -- but good enough
const std::vector<double> NPL_POS = {3984.91, -23.9022, 4963.31};
const std::vector<double> PTB_POS = {3843.9, 709.921, 5023.05}; // km
const std::vector<double> SYRTE_POS = {4202.72, 171.361, 4778.55};

const std::vector<std::vector<double>> LABPOS = {
    NPL_POS,   NPL_POS,   PTB_POS,   PTB_POS,  SYRTE_POS,
    SYRTE_POS, SYRTE_POS, SYRTE_POS, SYRTE_POS};

inline int GET_CLOCK_INDEX(std::string in_clock) {
  for (std::size_t i = 0; i < CLOCKS.size(); i++)
    if (in_clock == CLOCKS[i])
      return (int)i;
  std::cout << "\nFAILURE 71 in clockInfo: cant find " << in_clock << "\n\n";
  return -1; // bad error
}

} // namespace ClockInfo
