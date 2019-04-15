#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

bool readEachFile(std::string infn, std::vector<double> &T_list,
                  std::vector<double> &tint_list,
                  std::vector<std::vector<double>> &dX) {
  std::ifstream file;
  file.open(infn);
  {
    std::string line;
    // first line:
    std::getline(file, line);
    std::stringstream ssin(line);
    std::string junk;
    ssin >> junk; // first entry is junk
    while (true) {
      double tmp_T;
      ssin >> tmp_T;
      if (!ssin.good())
        break;
      T_list.push_back(tmp_T);
    }
    // Read subsequent lines: (dX, and tau_int)
    while (std::getline(file, line)) {
      std::stringstream ssin(line);
      double tmp_tint;
      ssin >> tmp_tint;
      tint_list.push_back(tmp_tint);
      std::vector<double> dX_tint;
      while (true) {
        double tmp_dX;
        ssin >> tmp_dX;
        if (!ssin.good())
          break;
        if (tmp_dX == 0)
          tmp_dX = 1. / 0.; // OK? ??
        dX_tint.push_back(tmp_dX);
      }
      if (dX_tint.size() != T_list.size()) {
        std::cout << "\nFAIL:61\n";
        return false;
      }
      dX.push_back(dX_tint);
    }
    if (dX.size() != tint_list.size()) {
      std::cout << "\nFAIL:65\n";
      return false;
    }
  }
  file.close();
  return true;
}

//******************************************************************************
int main(void) {

  // std::string infn = "coinc_dX_test1.txt";

  std::string label = "xx";

  std::vector<std::string> file_list; // = {"coinc_dX_test1.txt"};

  std::vector<double> T_list;
  std::vector<double> tint_list;
  std::vector<std::vector<double>> dX;

  std::vector<double> tint_min; // = {15.};
  std::vector<double> tint_max; // = {0};

  std::vector<double> T_out;    // = {1.,10.,100.};
  std::vector<double> tint_out; // = {15.,100., 1000.};

  // Read in the input file (lines marked with '!' are ignored)
  {
    std::ifstream ifs;
    ifs.open("combinedX.in");
    std::string line;
    while (std::getline(ifs, line)) {
      if (line.front() == '!')
        continue;
      label = line;
      break;
    }

    while (std::getline(ifs, line)) {
      if (line.front() == '!')
        continue;
      std::stringstream ssin(line);
      while (ssin.good()) {
        double tT;
        ssin >> tT;
        T_out.push_back(tT);
      }
      break;
    }

    while (std::getline(ifs, line)) {
      if (line.front() == '!')
        continue;
      std::stringstream ssin(line);
      while (ssin.good()) {
        double ttau;
        ssin >> ttau;
        tint_out.push_back(ttau);
      }
      break;
    }

    while (std::getline(ifs, line)) {
      if (line.front() == '!')
        continue;
      std::string tname;
      double tmin, tmax;
      std::stringstream ssin(line);
      ssin >> tname >> tmin >> tmax;
      file_list.push_back(tname);
      tint_min.push_back(tmin);
      tint_max.push_back(tmax);
    }
    ifs.close();
  }

  // return 1;
  for (size_t i = 0; i < file_list.size(); i++) {
    std::cout << file_list[i] << "\n";
  }

  // Loop through each file, constuct limit matrix:
  for (size_t i = 0; i < file_list.size(); i++) {
    std::string infn = file_list[i];
    std::vector<double> T_list_i;
    std::vector<double> tint_list_i;
    std::vector<std::vector<double>> dX_i;
    readEachFile(infn, T_list_i, tint_list_i, dX_i);
    if (tint_max[i] == 0)
      tint_max[i] = 1.e15; // XXX
    double tau_min = tint_min[i];
    double tau_max = tint_max[i];
    // std::cout<<tau_max<<"\n";
    for (size_t j = 0; j < dX_i.size(); j++) {
      double tint = tint_list_i[j];
      if (tint >= tau_min && tint <= tau_max)
        continue;
      for (size_t k = 0; k < dX_i[j].size(); k++) {
        dX_i[j][k] = 1. / 0.;
      }
    }
    if (dX.size() == 0) {
      dX = dX_i;
      tint_list = tint_list_i;
      T_list = T_list_i;
    } else {
      // make sure new files exactly the same as other:
      if (tint_list_i.size() != tint_list.size() ||
          tint_list_i.size() != tint_list.size() ||
          tint_list_i.front() != tint_list.front() ||
          tint_list_i.back() != tint_list.back() ||
          T_list_i.size() != T_list.size() ||
          T_list_i.front() != T_list.front() ||
          T_list_i.back() != T_list.back())
        std::cout << "\nFAIL:81 - " << infn << "\n";
      // if better, replace [but only if applicable tau_int region!]
      for (size_t j = 0; j < dX.size(); j++) {
        // double tint = tint_list[j];
        // if(tint < tau_min || tint > tau_max) continue;
        for (size_t k = 0; k < dX[0].size(); k++) {
          if (dX_i[j][k] > 1.)
            continue;
          if (dX_i[j][k] < dX[j][k] || dX[j][k] > 1.) {
            dX[j][k] = dX_i[j][k];
          }
        }
      }
    }
  }

  // find largest T for which we have timit:
  double Tmax = 0;
  for (int i = (int)T_list.size() - 1; i >= 0; --i) {
    Tmax = T_list[i];
    bool found = false;
    for (size_t j = 0; j < tint_list.size(); j++) {
      if (dX[j][i] < 1) {
        found = true;
        break;
      }
    }
    if (found)
      break;
  }
  // add to list:
  T_out.push_back(Tmax);

  double hcrv2 = 7.1e-6; // hbar*c*rho_DM*v_g^2 = in: (TeV/s)^2
  double hr_s = 60. * 60;

  // For each input T, make plot as fn of tau_int
  // Find T index: note: T_calc < T_target
  std::vector<int> iT_list;
  for (size_t i = 0; i < T_out.size(); i++) {
    double T_targ = T_out[i];
    int iT = 0;
    for (size_t i2 = 0; i2 < T_list.size(); i2++) {
      if (T_list[i2] >= T_targ - 0.01) {
        iT = int(i2);
        break;
      }
    }
    iT_list.push_back(iT);
  }

  // Output file: function of tau_int, for few T
  std::string of1name = "dX_tauint_" + label + ".txt";
  std::string ofLname = "Lambda_tauint_" + label + ".txt";
  std::ofstream of, of2;
  of.open(of1name);
  of2.open(ofLname);
  of.precision(2);
  of2.precision(2);
  of << "tau_int\\T ";
  of2 << "tau_int\\T ";
  for (size_t i = 0; i < iT_list.size(); i++) {
    of << T_list[iT_list[i]] << " ";
    of2 << T_list[iT_list[i]] << " ";
  }
  of << "\n";
  of2 << "\n";
  of.precision(4);
  of2.precision(4);
  for (size_t j = 0; j < tint_list.size(); j++) {
    of << tint_list[j] << " ";
    of2 << tint_list[j] << " ";
    for (size_t i = 0; i < iT_list.size(); i++) {
      int iT = iT_list[i];
      of << dX[j][iT] << " ";
      double Ts = T_list[iT] * hr_s; // convert to seconds:
      double tint = tint_list[j];
      of2 << sqrt(hcrv2 * Ts * tint / dX[j][iT]) << " ";
    }
    of << "\n";
    of2 << "\n";
  }
  of.close();
  of2.close();

  // For each input tau_int, make plot as fn of T
  // For this, don't bother with dX, just Lambda
  std::vector<int> itint_list;
  for (size_t i = 0; i < tint_out.size(); i++) {
    double tint_targ = tint_out[i];
    int itint = 0;
    for (size_t i2 = 0; i2 < tint_list.size(); i2++) {
      if (tint_list[i2] >= tint_targ - 0.1) {
        itint = int(i2);
        break;
      }
    }
    itint_list.push_back(itint);
  }

  std::string of2name = "Lambda_T_" + label + ".txt";
  of.open(of2name);
  of.precision(2);
  of << "T\\tau_int ";
  for (size_t i = 0; i < itint_list.size(); i++)
    of << tint_list[itint_list[i]] << " ";
  of << "\n";
  of.precision(4);
  for (size_t j = 0; j < T_list.size(); j++) {
    of << T_list[j] << " ";
    for (size_t i = 0; i < itint_list.size(); i++) {
      int itint = itint_list[i];
      double Ts = T_list[j] * hr_s; // convert to seconds:
      double tint = tint_list[itint];
      of << sqrt(hcrv2 * Ts * tint / dX[itint][j]) << " ";
    }
    of << "\n";
  }
  of.close();

  return 0;
}
