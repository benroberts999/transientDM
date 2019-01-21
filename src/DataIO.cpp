#include "DataIO.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace DataIO{

//******************************************************************************
template <typename T, typename U>
int read_text_XY(const std::string &in_fname, std::vector<T> &x,
    std::vector<U> &y, bool append)
/*
return 0 = all good.
return 1 = file is empty or doesn't exist (or can't be opened)
return 2 = something worse happened.. bug in program, or data file is messed up
*/
{
  std::ifstream ifs;
  ifs.open(in_fname.c_str());
  if(!ifs.good()) return 1;
  std::string str_line;
  if(!append){
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
  while(getline(ifs,str_line)){
    if(str_line.size()==0) continue;
    std::stringstream ss(str_line);
    T tmp_x=0;
    U tmp_y=0;
    ss>>tmp_x>>tmp_y;
    x.push_back(tmp_x);
    y.push_back(tmp_y);
  }
  ifs.close();
  if(x.size()!=y.size()) return 2;
  if(x.size()==0) return 1;
  return 0;
}

//Forward-declare templates:
template int read_text_XY(const std::string &in_fname, std::vector<double> &x,
   std::vector<double> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<float> &x,
   std::vector<float> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<int> &x,
   std::vector<int> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<int> &x,
   std::vector<double> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<int> &x,
   std::vector<float> &y, bool append);

}//namespace
