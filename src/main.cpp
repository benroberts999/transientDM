#include "ClockNetwork.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
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

  int tau_avg = 60;
  ClockNetwork net(filenames,tau_avg);

  long j_init = net._initial_epoch[3];
  long j_fin = j_init+6000;

  //make a net.printInfo() command!

  //return 1;
  int teff_min = 60; //must be greater than (or =) tau_avg! XXX
  int teff_max = 3000;
  int nteff = 9;
  //Define the tau_eff grid (exponentially spaced grid):
  int teff_ot0_min = teff_min/tau_avg;
  int teff_ot0_max = teff_max/tau_avg;
  if(nteff==1) teff_ot0_max = teff_ot0_min;
  std::vector<int> teff_ot0_list(nteff);
  teff_ot0_list[0] = teff_ot0_min;
  for(int i=1; i<nteff; i++){
    double x = double(i)/(nteff-1);
    double teff_tmp = teff_ot0_min*pow(double(teff_ot0_max)/teff_ot0_min,x);
    int iteff_tmp = (int)round(teff_tmp);
    //this is to stop "doubles" from occuring (teff is an integer):
    if(i>0 && iteff_tmp<=teff_ot0_list[i-1]) iteff_tmp = teff_ot0_list[i-1]+1;
    teff_ot0_list[i] = iteff_tmp;
  }

  //for(int i=1; i<15; i++){
  for(size_t it=0; it<teff_ot0_list.size(); it++){
    int tau_int_on_tau = teff_ot0_list[it];
    std::cout<<"tau_eff="<<tau_int_on_tau*tau_avg<<"\n";

    std::vector<std::vector<double> > s;
    int window_tauint_mult = 1;
    net.genSignalTemplate(s,window_tauint_mult,tau_int_on_tau);

    int j_step = tau_avg; //XXX How to determine Tobs?
    int num_j0 = 0;
    for(long jbeg = j_init; jbeg<=j_fin; jbeg+=j_step){
      // std::cerr<<jbeg<<"\n";

      int Jw = s[0].size();
      int max_bad = 0;

      std::vector<int> indep_pairs;
      net.formIndependentSubnet(indep_pairs,jbeg,Jw,max_bad);

      int N_pairs = (int)indep_pairs.size();
      //std::cout<<"Idep pairs: "<<N_pairs<<"\n";
      if(N_pairs<1) continue;

      ++num_j0;

      auto dHs_sHs = net.calculate_dHs_sHs(indep_pairs,s,jbeg);
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
