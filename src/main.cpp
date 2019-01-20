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

  std::string
  fn="./data/cppdiff_SYRTEHg-PTBSr.dat";
  net.readInDataFile(fn,1);
  fn="./data/cppdiff_PTBSr-PTBYb.dat";
  net.readInDataFile(fn,1);

  //std::cout<<"X: "<<net.new_vector3[0]<<"\n";

  std::cout<<net._K_AB[0]<<"\n";
  std::cout<<net._delta_omega[0][1]<<"\n";
  std::cout<<net._initial_time[0]<<"\n";

  std::cout<<net._K_AB[1]<<"\n";
  std::cout<<net._delta_omega[1][1]<<"\n";
  std::cout<<net._initial_time[1]<<"\n";

  net.calculateSigma0();
  std::cout<<net._sigma0[0]<<" "<<net._sigma0[1]<<"\n";

  net.rankClockPairs();

  std::vector<std::vector<double> > s;
  net.genSignalTemplate(s,1,2);
  for(size_t i=0; i<s.size(); i++){
    for(size_t j=0; j<s[i].size(); j++){
      std::cout<<s[i][j]/net._K_AB[i]<<" ";
    }
    std::cout<<"\n";
  }

  std::vector<int> indep_pairs;
  net.formIndependentSubnet(indep_pairs,566818,10,0);
  for(auto i : indep_pairs)
    std::cout<<i<<"\n";

}
