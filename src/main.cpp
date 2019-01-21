#include "ClockNetwork.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "DataIO.h"
#include "ChronoTimer.h"
#include <array>


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
  ChronoTimer timer;
  timer.start();

  std::vector<std::string> filenames;
  getFileNames(filenames,"./clocklist.in");

  int tau_avg = 60;
  int min_N_pairs = 2; //XXX input! Must be greater than 0
  int max_bad_avg = 0;
  int iday = 0; //in days
  int fday = 50;
  int teff_min = 60; //must be greater than (or =) tau_avg! XXX
  int teff_max = 10000;
  int nteff = 32;
  int nJeffW = 1; //window multiplier
  bool force_PTB_SrYb=false;
  bool force_SyrHbNplYb=false;

  timer.start();
  ClockNetwork net(filenames,tau_avg,max_bad_avg);
  std::cout<<"Time to read files: "<<timer.lap_reading_str()<<"\n\n";

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


  //Define the tau_int = tau_eff grint to search:
  //XXX Define prperly (equation) link between tau_int and tau_eff !!!
  //note: j_eff := tau_eff / tau_0 = tau_eff / tau_avg
  int jeff_min = teff_min/tau_avg;
  int jeff_max = teff_max/tau_avg;

  //Define the tau_eff grid (exponentially spaced grid):
  std::vector<int> jeff_grid;
  defineIntegerLogGrid(jeff_grid, jeff_min, jeff_max, nteff);

  std::vector<double> da_max(nteff);
  std::vector<double> Del_da(nteff);
  std::vector<double> R_max(nteff);
  std::vector<double> T_obs(nteff);

  // #pragma omp parallel for
  for(size_t it=0; it<jeff_grid.size(); it++){
    int j_eff = jeff_grid[it];

    std::vector<std::vector<double> > s;
    s.reserve(N_tot_pairs);

    // Jw = nJeffW*j_eff, or nJeffW*j_eff+1  (Jw must be odd)
    net.genSignalTemplate(s,nJeffW,j_eff, TDProfile::Gaussian);
    int Jw = (int)s[0].size();

    int j_step = 1;
    //OR j_eff, or j_eff/2 (minimum of 1, must be int!)
    //j_step = j_eff;  //XXX make an input option!! XXX
    //if(j_step==0) j_step = 1;

    int num_j_used = 0;
    for(long jbeg = j_init; jbeg<=j_fin; jbeg+=j_step){

      int max_bad_inJw = 0;

      std::vector<int> indep_pairs;
      net.formIndependentSubnet(indep_pairs,jbeg,Jw,max_bad_inJw,
        force_PTB_SrYb,force_SyrHbNplYb);
      //XXX Update this to a "must include" list ?

      int N_pairs = (int)indep_pairs.size();
      if(N_pairs<min_N_pairs) continue;

      ++num_j_used;

      auto xHs = net.calculate_dHs_sHs(indep_pairs,s,jbeg);

      double da_bf = fabs(xHs.dHs/xHs.sHs);

      if(da_bf>da_max[it]){
        da_max[it] = da_bf;
        Del_da[it] = 1./sqrt(xHs.sHs);
      }

    }

    double Tobs_hr = num_j_used*j_step*tau_avg/(60*60.);
    T_obs[it] = Tobs_hr;

    std::cout<<j_eff*tau_avg<<" "<<(da_max[it] + Del_da[it])
      <<"   (Tobs = "<<T_obs[it]<<" hr)"<<"\n";

  }

  std::cout<<"\nTotal time: "<<timer.reading_str()<<"\n";
  return 0;
}//End main()
