#include "ClockNetwork.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "DataIO.h"



void getFileNames(std::vector<std::string> &fnames,
  const std::string &input_fn)
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


int main(){

  std::vector<std::string> filenames;
  getFileNames(filenames,"./clocklist.in");

  int tau_avg = 500;
  ClockNetwork net(filenames,tau_avg);

  //make a net.printInfo() command!

  //return 1;

  for(int i=1; i<15; i++){
    int tau_int_on_tau = i;

    std::vector<std::vector<double> > s;
    int window_tauint_mult = 1;
    net.genSignalTemplate(s,window_tauint_mult,tau_int_on_tau);

    long j_beg_epoch_0 = net._initial_epoch[3];
    long j_beg_epoch=j_beg_epoch_0;
    //for(long j_beg_epoch=j_beg_epoch_0; j_beg_epoch<j_beg_epoch_0+1000; j_beg_epoch++)
    {
      int Jw = s[0].size();
      int max_bad = 0;

      std::vector<int> indep_pairs;
      net.formIndependentSubnet(indep_pairs,j_beg_epoch,Jw,max_bad);

      int N_pairs = (int)indep_pairs.size();
      //std::cout<<"Idep pairs: "<<N_pairs<<"\n";
      if(N_pairs<1) continue;

      auto dHs_sHs = net.calculate_dHs_sHs(indep_pairs,s,j_beg_epoch);
      std::cout<<dHs_sHs[0]/dHs_sHs[1]<<"\n";
    }

  }


  //
  //
  //
  //
  // std::cout<<net._sigma0[0]<<" "<<net._sigma0[1]<<"\n";
  //
  //
  //
  // std::vector<std::vector<double> > s;
  // net.genSignalTemplate(s,1,3);
  // for(size_t i=0; i<s.size(); i++){
  //   for(size_t j=0; j<s[i].size(); j++){
  //     std::cout<<s[i][j]/net._K_AB[i]<<" ";
  //   }
  //   std::cout<<"\n";
  // }
  //
  // for(auto i : net._ranked_index_list)
  //   std::cout<<net._initial_epoch[i]<<"\n";
  //
  // std::vector<int> indep_pairs;
  // net.formIndependentSubnet(indep_pairs,26829,10,0);
  // for(auto i : indep_pairs){
  //   std::cout<<i<<"\n";
  //   std::cout<<net._clock_name_A[i]<<"-"<<net._clock_name_B[i]<<"\n";
  // }
  //
  // auto dHs_sHs = net.calculate_dHs_sHs(indep_pairs,s,26829);
  // std::cout<<dHs_sHs[0]/dHs_sHs[1]<<"\n";
  //
  // for(int i=1; i<15; i++){
  //   int tau_int_on_tau = i;
  //
  //   std::vector<std::vector<double> > s;
  //   net.genSignalTemplate(s,2,tau_int_on_tau);
  //
  //   std::vector<int> indep_pairs;
  //   net.formIndependentSubnet(indep_pairs,26829,s[0].size(),0);
  //
  //   std::cout<<"Idep pairs: "<<indep_pairs.size()<<"\n";
  //
  //   auto dHs_sHs = net.calculate_dHs_sHs(indep_pairs,s,26829);
  //   std::cout<<dHs_sHs[0]/dHs_sHs[1]<<"\n";
  //
  // }

}
