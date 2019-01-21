#pragma once
#include <vector>
#include <string>

namespace DataIO{

  template <typename T, typename U>
  int read_text_XY(const std::string &in_fname, std::vector<T> &x,
      std::vector<U> &y, bool append=false);

}
