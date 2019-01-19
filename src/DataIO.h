#pragma once
#include <vector>
#include <string>

namespace DataIO{

  template <typename T, typename U>
  int read_text_XY(const std::string &in_fname, std::vector<T> &x,
      std::vector<U> &y, bool append=false);

  // template <typename T>
  // int writeXYTextFile(const std::string &out_fname,
  //   const std::vector<double> &x, const std::vector<double> &y,
  //   bool append=false, int precision=10);
  //
  // template <typename T>
  // int writeYTextFile(const std::string &out_fname, const std::vector<T> &y,
  //   bool append=false, int precision=10);

}
