#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace DataIO{

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
  std::string path="";
  std::string str;
  bool header = true;
  while(getline(ifs,str)){
    if(str.substr(0,1)=="!"||str.substr(0,1)=="#") continue;
    if(header){
      path = str;
      header = false;
    }else{
      fnames.push_back(path+str);
    }
  }
  ifs.close();
}

//******************************************************************************
template <typename T, typename U>
int read_text_XY(const std::string &in_fname, std::vector<T> &x,
      std::vector<U> &y, bool append=false)
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

}
