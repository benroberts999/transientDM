#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace DataIO {

//******************************************************************************
inline void getFileNames(std::vector<std::string> &fnames,
                         const std::string &input_fn)
/*
Opens an input text files (called: input_fn)
Reads each line into an array.
First line of text-file is an absolute (or relative) path name.
Each subsequent line contains just the file names.
Any line that begines with '!' or '#' is ignored
*/
{
  std::ifstream ifs;
  ifs.open(input_fn);
  std::string path = "";
  std::string str;
  bool header = true;
  while (getline(ifs, str)) {
    if (str.substr(0, 1) == "!" || str.substr(0, 1) == "#")
      continue;
    if (header) {
      path = str;
      header = false;
    } else {
      fnames.push_back(path + str);
    }
  }
  ifs.close();
}

//******************************************************************************
template <typename T, typename U>
int read_text_XY(const std::string &in_fname, std::vector<T> &x,
                 std::vector<U> &y, bool append = false)
/*
return 0 = all good.
return 1 = file is empty or doesn't exist (or can't be opened)
return 2 = something worse happened.. bug in program, or data file is messed up
*/
{
  std::ifstream ifs;
  ifs.open(in_fname.c_str());
  if (!ifs.good())
    return 1;
  std::string str_line;
  if (!append) {
    x.clear();
    y.clear();
  }
  // //********************************** ???
  // long number_of_lines = 0;
  // while (std::getline(ifs, str_line))
  //       ++number_of_lines;
  // ifs.clear();
  // ifs.seekg(0, std::ios::beg);
  // //////////
  // x.reserve(x.size()+number_of_lines);
  // y.reserve(y.size()+number_of_lines);
  // //**********************************
  while (getline(ifs, str_line)) {
    if (str_line.size() == 0)
      continue;
    std::stringstream ss(str_line);
    T tmp_x = 0;
    U tmp_y = 0;
    ss >> tmp_x >> tmp_y;
    x.push_back(tmp_x);
    y.push_back(tmp_y);
  }
  ifs.close();
  if (x.size() != y.size())
    return 2;
  if (x.size() == 0)
    return 1;
  return 0;
}

//******************************************************************************
template <typename T, typename U>
int write_text_XY(const std::string &out_fname, std::vector<T> &x,
                  std::vector<U> &y, bool append = false)
/*
return 0 = all good.
return 1 = file is empty or doesn't exist (or can't be opened)
return 2 = something worse happened.. bug in program, or data file is messed up
*/
{
  auto wr_mode = append ? std::ios_base::app : std::ios_base::out;
  std::ofstream ofs(out_fname.c_str(), wr_mode);

  if (x.size() != y.size()) {
    std::cerr << "\nFail 96 write_text_XY: vectors unequal length!\n";
    return 2;
  }

  for (size_t i = 0; i < x.size(); i++)
    ofs << x[i] << " " << y[i] << "\n";

  return 0;
}

//******************************************************************************
/*
Uses compile-time recursion to get access to elements of tuple.
Specifically, string-streams data from a string vector into tuple.
Works with a tuple of references, i.e., std::forward_as_tuple
Idea from:
https://stackoverflow.com/questions/1198260/iterate-over-tuple/23142715
*/
template <std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
stringstreamVectorIntoTuple(std::vector<std::string>, std::tuple<Tp...> &) {}

template <std::size_t I = 0, typename... Tp>
    inline typename std::enable_if <
    I<sizeof...(Tp), void>::type
    stringstreamVectorIntoTuple(std::vector<std::string> lst,
                                std::tuple<Tp...> &t) {
  if (I > lst.size())
    std::cerr << "\nFAIL 34 in FileIO: list shorter than tuple\n";
  std::stringstream(lst[I]) >> std::get<I>(t);
  stringstreamVectorIntoTuple<I + 1, Tp...>(lst, t);
}

//******************************************************************************
inline std::vector<std::string> readInputFile(const std::string &fname) {
  std::vector<std::string> entry_list;
  std::ifstream file(fname);
  std::string line;
  while (getline(file, line)) {
    std::stringstream ss(line);
    std::string entry;
    while (ss >> entry) {
      if (entry.at(0) == '!' || entry.at(0) == '#')
        break;
      else
        entry_list.push_back(entry);
    }
  }
  return entry_list;
}

//******************************************************************************
template <typename... Tp>
void setInputParameters(std::string infile, std::tuple<Tp...> &tp) {
  auto input = readInputFile(infile);
  stringstreamVectorIntoTuple(input, tp);
}

} // namespace DataIO
