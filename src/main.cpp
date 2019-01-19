#include "ClockNetwork.h"
#include <iostream>
#include <vector>
#include <string>
#include "DataIO.h"


int main(){

  std::cout<<"hello\n";

  //ClockNetwork net();
  std::vector<double> x,y; //{1,2,3,4,5,6,7,8,9,10};
  // std::vector<double> y
  DataIO::read_text_XY("input.txt",x,y);

  for(auto xi : x)
    std::cout<<xi<<"\n";

  ClockNetwork net;

  std::string fn="./data/cppdiff_PTBSr-PTBYb.dat";
  net.readInDataFile(fn);

  std::cout<<net._K_AB[0]<<"\n";

  std::cout<<net._delta_omega[0][1]<<"\n";
  std::cout<<net._initial_time[0]<<"\n";

}
