#include "DataIO.h"
// #include <cmath>
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
  // x.reserve(50000);
  // y.reserve(50000);
  std::ifstream ifs;
  ifs.open(in_fname.c_str());
  if(!ifs.good()) return 1;
  std::string str_line;
  if(!append){
    x.clear();
    y.clear();
  }
  /////
  // long number_of_lines = 0;
  // while (std::getline(ifs, str_line))
  //       ++number_of_lines;
  // ifs.clear();
  // ifs.seekg(0, std::ios::beg);
  // //////////
  // x.reserve(x.size()+number_of_lines);
  // y.reserve(y.size()+number_of_lines);
  // //////////
  while(getline(ifs,str_line)){
    if(str_line.size()==0) continue;
    std::stringstream ss(str_line);
    T tmp_x=0;
    U tmp_y=0;
    ss>>tmp_x>>tmp_y;
    x.push_back(tmp_x);
    y.push_back(tmp_y);
    // x.emplace_back(tmp_x);
    // y.emplace_back(tmp_y);
  }
  ifs.close();
  if(x.size()!=y.size()) return 2;
  if(x.size()==0) return 1;
  return 0;
}
template int read_text_XY(const std::string &in_fname, std::vector<double> &x, std::vector<double> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<float> &x, std::vector<float> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<int> &x, std::vector<int> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<int> &x, std::vector<double> &y, bool append);
template int read_text_XY(const std::string &in_fname, std::vector<int> &x, std::vector<float> &y, bool append);

// //******************************************************************************
// int writeXYTextFile(const std::string &out_fname, const std::vector<double> &x,
//     const std::vector<double> &y, bool append, int precision)
// /*
// precision = 10 by default.
// return 0 = all good.
// return 1 = no data to write out
// return 2 = Bad. Cannot create file, or x,y of different length.
// */
// {
//   if(x.size()!=y.size()) return 2;
//   if(x.size()==0) return 1;
//   std::ofstream ofs;
//   if(append) ofs.open(out_fname, std::ios_base::app);
//   else ofs.open(out_fname.c_str());
//   if(!ofs.good()) return 2;
//   ofs.precision(precision);
//   //if(y.size()==0) ofs<<"\n";
//   for(size_t i=0; i<x.size(); i++){
//     ofs<<x[i]<<" "<<y[i]<<"\n";
//   }
//   ofs.close();
//   return 0;
// }
//
// //******************************************************************************
// int writeYTextFile(const std::string &out_fname, const std::vector<double> &y,
//     bool append, int precision)
// /*
// precision = 10 by default.
// return 0 = all good.
// return 1 = no data to write out
// return 2 = Bad. Cannot create file, or x,y of different length.
// */
// {
//   if(y.size()==0) return 1;
//   std::ofstream ofs;
//   if(append) ofs.open(out_fname, std::ios_base::app);
//   else ofs.open(out_fname.c_str());
//   if(!ofs.good()) return 2;
//   ofs.precision(precision);
//   //if(y.size()==0) ofs<<"\n";
//   for(size_t i=0; i<y.size(); i++){
//     ofs<<y[i]<<"\n";
//   }
//   ofs.close();
//   return 0;
// }

}//namespace
