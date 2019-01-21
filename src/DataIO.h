#pragma once
#include <vector>
#include <string>

namespace DataIO{

  void getFileNames(std::vector<std::string> &fnames,
    const std::string &input_fn);

  template <typename T, typename U>
  int read_text_XY(const std::string &in_fname, std::vector<T> &x,
      std::vector<U> &y, bool append=false);

}
