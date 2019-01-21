#include "ClockNetwork.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "DataIO.h"
#include "ChronoTimer.h"


//******************************************************************************
void getFileNames(std::vector<std::string> &fnames,
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
void defineIntegerLogGrid(std::vector<int> &grid, int min, int max, int N)
/*
Forms a logarithmically-spaced grid of integers, between [min,max] in N steps.
Note: Because each point is an integer, if the step-size is too small,
you run the risk of getting 'doubles' of some grid points (same point more than
once).
I fix this, by just incrementing by 1 in these cases.
Technically, this means it's possible for the actual last grid point to go
beyond the given 'max' (if N > (max-min)), and the grid will then just be linear
But, it _is_ guarenteed to be exactly N elements long
*/
{
  grid.resize(N);
  grid[0] = min;

  for(int i=1; i<N; i++){
    double x = double(i)/(N-1);
    double y = min*pow(double(max)/min,x);
    int z = (int) round(y);
    //this is to stop "doubling up" from occuring:
    if(i>0 && z<=grid[i-1]) z = grid[i-1]+1;
    grid[i] = z;
  }

}



//******************************************************************************
int main(){

  std::vector<std::string> filenames;
  getFileNames(filenames,"./clocklist.in");

  int tau_avg = 60;
  int min_N_pairs = 2; //XXX input! Must be greater than 0
  int max_bad_avg = 0;
  int iday = 0; //in days
  int fday = 50;
  int teff_min = 60; //must be greater than (or =) tau_avg! XXX
  int teff_max = 5000;
  int nteff = 32;
  int nJeffW = 1; //window multiplier



  ChronoTimer timer;
  timer.start();
  ClockNetwork net(filenames,tau_avg,max_bad_avg);
  std::cout<<"Timer: "<<timer.lap_reading_str()<<"\n\n";

  int N_tot_pairs = net.get_NtotPairs();
  if(N_tot_pairs < min_N_pairs){
    std::cerr<<"No clocks? Check data file paths\n";
    return 1;
  }

  //initial/final EPOCHS to search.
  //Epochs are time/tau_0.
  //Initial/final: convert days -> seconds
  long j_init = (long) (iday*24*60*60 / tau_avg); //in epochs!
  long j_fin = (long) ((fday*24*60*60 + 1) / tau_avg);
  //
  std::cout<<"Running for days (since MJD:"<<MJD_DAY_ZERO<<"): "
    <<j_init*tau_avg/(24*60*60)<<" -> "<<j_fin*tau_avg/(24*60*60)<<"\n";
  std::cout<<"= "<<j_init*tau_avg<<" -> "<<j_fin*tau_avg<<" s\n";
  std::cout<<"= "<<j_init<<" -> "<<j_fin<<" epochs (w/ tau_0=tau_avg = "
    <<tau_avg<<")\n\n";
  // return 1;


  //Define the tau_int = tau_eff grint to search:
  //XXX Define prperly (equation) link between tau_int and tau_eff !!!
  //note: j_eff := tau_eff / tau_0 = tau_eff / tau_avg
  int jeff_min = teff_min/tau_avg;
  int jeff_max = teff_max/tau_avg;

  //Define the tau_eff grid (exponentially spaced grid):
  std::vector<int> jeff_grid;
  defineIntegerLogGrid(jeff_grid, jeff_min, jeff_max, nteff);

  // #pragma omp parallel for
  for(size_t it=0; it<jeff_grid.size(); it++){
    int j_eff = jeff_grid[it];

    std::vector<std::vector<double> > s;
    s.reserve(N_tot_pairs);

    // Jw = nJeffW*j_eff, or nJeffW*j_eff+1  (Jw must be odd)
    //XXX make above an input option!! XXX
    net.genSignalTemplate(s,nJeffW,j_eff, TDProfile::Gaussian);
    int Jw = s[0].size();

    double max_da = 0;
    double Del_da = 0;

    int j_step = 1;//tau_avg; XXX
    //OR j_eff, or j_eff/2 (minimum of 1, must be int!)
    j_step = j_eff/2;
    if(j_step==0) j_step = 1;

    int num_j_used = 0;
    for(long jbeg = j_init; jbeg<=j_fin; jbeg+=j_step){

      int max_bad = 0;

      std::vector<int> indep_pairs;
      net.formIndependentSubnet(indep_pairs,jbeg,Jw,max_bad);

      int N_pairs = (int)indep_pairs.size();
      if(N_pairs<min_N_pairs) continue;

      ++num_j_used;

      auto xHs = net.calculate_dHs_sHs(indep_pairs,s,jbeg);
      // std::cout<<it<<" "<<jbeg<<" "<<N_pairs<<"\n";
      // std::cout<<" --> "<<xHs.dHs/xHs.sHs<<"\n";

      double da_bf = fabs(xHs.dHs/xHs.sHs);

      if(da_bf>max_da){
        max_da = da_bf;
        Del_da = 1./sqrt(xHs.sHs);
      }

    }

    std::cout<<j_eff*tau_avg<<" "<<(max_da + Del_da)
     <<"   (Tobs = "<<num_j_used*j_step*tau_avg/(60*60.)<<" hr)"<<"\n";

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
