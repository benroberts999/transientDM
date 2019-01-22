#include "ClockNetwork.h"
#include "DataIO.h"
#include "ChronoTimer.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

enum class WhichOutput{ limit, R,  limit_and_R, Rthresh};

void defineIntegerLogGrid(std::vector<int> &grid, int min, int max, int N);
void defineDoubleLogGrid(std::vector<double> &grid, double min, double max,
  int N);

void dmSearch_tau_int(
  const ClockNetwork &net, int teff_min, int teff_max, int nteff,
  TDProfile profile, int nJeffW, double n_sig, double iday, double fday,
  int min_N_pairs, bool force_PTB_SrYb, bool force_SyrHbNplYb,
  std::string olabel, WhichOutput whichoutput);

void outputConstraints( const std::vector<int> &jeff_grid, int tau_avg,
  const std::vector<double> &da_max, const std::vector<double> &Del_da,
  const std::vector<double> &T_obs, double n_sig, const std::string &olabel);

void outputR_teff(const std::vector<int> &jeff_grid, int tau_avg,
  const std::vector<double> &R_max,
  const std::vector<double> &da_max, const std::vector<double> &Del_da,
  const std::vector<double> &T_obs, const std::string &olabel);


//******************************************************************************
int main(){
  ChronoTimer timer;

  //Input parameters (to be read from file)
  std::string clock_list_infn; //File that contains clock list
  double iday, fday; //First/last day
  int tau_avg, max_bad_avg; //average data?
  int teff_min, teff_max, nteff;  //tau_eff grid
  int min_N_pairs;  //min. allowable # of clock pairs
  bool force_PTB_SrYb, force_SyrHbNplYb;  //Force specific clocks?
  TDProfile profile;  //TD profile (flat/gaussian)
  int n_sig; //[1=68%CL, 1.645=90%, 2=95.5%]
  std::string olabel;
  WhichOutput whichoutput; //enum for what output/calculation to do

  int nJeffW = 1; //window multiplier (not read from file)

  //Read in from input file:
  {
    std::ifstream ifs;
    std::string jnk;
    int iPSY, iSHNY, iP, iwhat;
    ifs.open("transientDM.in");
    ifs >> clock_list_infn;                 getline(ifs,jnk);
    ifs >> iday >> fday;                    getline(ifs,jnk);
    ifs >> tau_avg >> max_bad_avg;          getline(ifs,jnk);
    ifs >> teff_min >> teff_max >> nteff;   getline(ifs,jnk);
    ifs >> min_N_pairs >> iPSY >> iSHNY;    getline(ifs,jnk);
    ifs >> iP;                              getline(ifs,jnk);
    ifs >> n_sig;                           getline(ifs,jnk);
    ifs >> iwhat;                           getline(ifs,jnk);
    ifs >> olabel;                          getline(ifs,jnk);
    ifs.close();
    force_PTB_SrYb   = (iPSY ==1) ? true : false;
    force_SyrHbNplYb = (iSHNY==1) ? true : false;
    profile = (iP==1) ? TDProfile::Flat : TDProfile::Gaussian;
    olabel = (olabel=="na") ? "" : "-"+olabel;
    if(iwhat==0) whichoutput = WhichOutput::limit;
    else if(iwhat==1) whichoutput = WhichOutput::R;
    else if(iwhat==2) whichoutput = WhichOutput::limit_and_R;
    else if(iwhat==3) whichoutput = WhichOutput::Rthresh;
    else whichoutput = WhichOutput::limit; //default..
  }

  //Some checks for valid input data:
  if(teff_min<tau_avg){
    std::cerr<<"teff_min<tau_avg\n";
    return 1;
  }
  if(teff_max==teff_min) nteff=1;
  if(nteff==1) teff_max=teff_min;
  if(min_N_pairs<=0) min_N_pairs=1;

  //Read-in the list of clock-data filenames
  std::vector<std::string> filenames;
  DataIO::getFileNames(filenames,clock_list_infn);

  //read in the clock data + instantiate ClockNetwork object:
  timer.start();
  ClockNetwork net(filenames,tau_avg,max_bad_avg);
  std::cout<<"Time to read files: "<<timer.lap_reading_str()<<"\n\n";

  int N_tot_pairs = net.get_NtotPairs();
  if(N_tot_pairs < min_N_pairs){
    std::cerr<<"No clocks? Check data file paths\n";
    return 1;
  }

  //Output some info to screen:
  std::cout<<"Running for days (since MJD:"<<CNconsts::MJD_DAY_ZERO<<"): "
    <<iday<<" -> "<<fday<<" = "<<iday*24*60*60<<" -> "<<fday*24*60*60<<" s.\n";
  std::cout<<"With: tau_0 = tau_avg = " <<tau_avg<<"s.\n";
  std::cout<<"There are a total of: "<<N_tot_pairs<<" clock pairs. "
    <<" Only running analysis for N>"<<min_N_pairs<<"\n";
  std::cout<<"Searching between tau_int = "<<teff_min<<"->"<<teff_max
    <<"s, in: "<<nteff<<" steps.\n";
  if(profile==TDProfile::Gaussian)
    std::cout<<"Assuming a Gaussian profile.\n";
  else
    std::cout<<"Assuming a flat profile.\n";
  std::cout<<"\n";

  // std::cout<<"\nRandomising:\n"<<std::flush;
  // net.replaceWithRandomNoise(FillGaps::yes);
  // std::cout<<"Randomised.\n"<<std::flush;

  timer.start();
  dmSearch_tau_int(net, teff_min,teff_max,nteff, profile,nJeffW,n_sig,
    iday,fday, min_N_pairs,force_PTB_SrYb,force_SyrHbNplYb,olabel,
    whichoutput); //yay! So few parameters!
  std::cout<<"Time to analyse data: "<<timer.lap_reading_str()<<"\n";

  std::cout<<"\nTotal time: "<<timer.reading_str()<<"\n";
  return 0;
}//END main()
//******************************************************************************


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
void defineDoubleLogGrid(std::vector<double> &grid, double min, double max,
  int N)
/*
Forms a logarithmically-spaced grid of doubles, between [min,max] in N steps.
*/
{
  grid.clear();
  grid.reserve(N);
  grid.push_back(min);
  for(int i=1; i<N; i++){
    double x = double(i)/(N-1);
    grid.emplace_back(min*pow(max/min,x));
  }
}


//******************************************************************************
void dmSearch_tau_int(
  const ClockNetwork &net, int teff_min, int teff_max, int nteff,
  TDProfile profile, int nJeffW, double n_sig, double iday, double fday,
  int min_N_pairs, bool force_PTB_SrYb, bool force_SyrHbNplYb,
  std::string olabel, WhichOutput whichoutput)
/*

*/
{

  bool outputLimit = true;
  bool outputR = false;
  switch(whichoutput){
    case WhichOutput::limit :
      outputLimit = true;
      outputR = false;
      break;
    case WhichOutput::R :
      outputLimit = false;
      outputR = true;
      break;
    case WhichOutput::limit_and_R :
      outputLimit = true;
      outputR = true;
      break;
    case WhichOutput::Rthresh :
      outputLimit = false;
      outputR = false;
      std::cout<<"Rthresh Not implemented yet!\n";
      return;
  }

  int tau_avg = net.get_tau0();
  int max_bad_inJw = 0; //just leave as zero?

  //initial/final epochs to search. (Epochs are time/tau_0.)
  //Convert days -> seconds -> epochs
  long j_init = (long) (iday*24*60*60 / tau_avg); //in epochs!
  long j_fin = (long) ((fday*24*60*60) / tau_avg) + 1;

  //Convert tau_eff (=tau_int) to j_eff (epoch)
  int jeff_min = teff_min/tau_avg;
  int jeff_max = teff_max/tau_avg;
  //Define the tau_eff grid (exponentially spaced grid):
  std::vector<int> jeff_grid;
  defineIntegerLogGrid(jeff_grid, jeff_min, jeff_max, nteff);

  //Arrays to store results:
  std::vector<double> da_max(nteff);
  std::vector<double> Del_da(nteff);
  std::vector<double> R_max(nteff);
  std::vector<double> T_obs(nteff);

  int N_tot_pairs = net.get_NtotPairs();

  #pragma omp parallel for
  for(size_t it=0; it<jeff_grid.size(); it++){
    int j_eff = jeff_grid[it];

    //signal template:
    std::vector<std::vector<double> > s;
    s.reserve(N_tot_pairs);
    net.genSignalTemplate(s,nJeffW,j_eff, profile);
    int Jw = (int)s[0].size(); //store window

    //Loop over j0 (t0): note: jbeg is _beginning_ of Jw window!
    long j_step = 1; //or tau_eff? tau_eff/2? Makes little diff!
    long num_j_used = 0; //count
    for(long jbeg = j_init; jbeg<=j_fin; jbeg+=j_step){

      //Form the independent sub-network:
      std::vector<int> indep_pairs;
      net.formIndependentSubnet(indep_pairs,jbeg,Jw,max_bad_inJw,
        force_PTB_SrYb,force_SyrHbNplYb);
      int N_pairs = (int)indep_pairs.size();
      if(N_pairs<min_N_pairs) continue;

      ++num_j_used; //count 'good'/used data points (for Tobs)

      auto xHs = net.calculate_dHs_sHs(indep_pairs,s,jbeg);

      //Maximise da and R (sepperately!)
      double da_bf = fabs(xHs.dHs)/xHs.sHs;
      double Delta = 1./sqrt(xHs.sHs);
      double R     = fabs(xHs.dHs)/sqrt(2*xHs.sHs);
      if(da_bf+Delta>da_max[it]+Del_da[it])
      //if(da_bf>da_max[it])
      {
        da_max[it] = da_bf;
        Del_da[it] = Delta; //1./sqrt(xHs.sHs);
      }
      if(R>R_max[it])
        R_max[it] = R;

    }// j (t0) loop

    //store observation time (for this tau_eff)
    double Tobs_hr = double(num_j_used*j_step*tau_avg)/(60*60.);
    T_obs[it] = Tobs_hr;

  }// tau_eff

  // for tau_eff = 60,120,1020,10020
  // => j_eff = 1,2,17,167
  std::cout<<"Summary of results:\n";
  std::cout<<"tau_eff  da_max    Da         Rmax   Tobs\n";
  for(size_t i=0; i<jeff_grid.size(); i++){
    auto jf = jeff_grid[i];
    if(jf==1 || jf==2 || jf==17 || jf==167)
      printf("%5is   %7.1e   %7.1e  %6.2f   %5.2fhr\n",
      jf*tau_avg,da_max[i],Del_da[i],R_max[i],T_obs[i]);
  }
  std::cout<<"\n";

  if(outputLimit)
    outputConstraints(jeff_grid,tau_avg,da_max,Del_da,T_obs,n_sig,olabel);

  if(outputR)
    outputR_teff(jeff_grid,tau_avg,R_max,da_max, Del_da,T_obs, olabel);

}


//******************************************************************************
void outputConstraints(
  const std::vector<int> &jeff_grid, int tau_avg,
  const std::vector<double> &da_max,
  const std::vector<double> &Del_da,
  const std::vector<double> &T_obs,
  double n_sig, const std::string &olabel)
{
  //Define output grids for limits.
  //(hard-coded output grids!)
  int odim_T = 128;
  double oT_min = 0.5; //hours
  double oT_max = 50.; //hours
  std::vector<double> T_grid;
  defineDoubleLogGrid(T_grid,oT_min,oT_max,odim_T);
  int odim_tint = 2048;
  double otint_min = 1.; //seconds
  double otint_max = 1.e5; //seconds
  std::vector<double> tint_grid;
  defineDoubleLogGrid(tint_grid,otint_min,otint_max,odim_tint);

  //Vector to store output limits
  std::vector<std::vector<double> > dX(odim_T, std::vector<double>(odim_tint));

  //Poisson factor;
  // require: T_obs > T / pois_fac
  // (or, T_max = T_obs * pois_fac)
  double pois_fac = -1./log(1.-erf(0.7*n_sig));

  //Sort the limits from the tau_eff grid into the tau_int grid
  //Taking observation time into account for T
  //Note: tau_int << T (otherwise, not transient! take 2x?)
  for(size_t iT=0; iT<T_grid.size(); iT++){
    double T = T_grid[iT];
    for(size_t iteff=0; iteff<jeff_grid.size(); iteff++){
      double t_eff = jeff_grid[iteff]*tau_avg;
      double Tobs = T_obs[iteff];
      if(T > Tobs*pois_fac) continue; //could probs beak, but meh
      double da = da_max[iteff] + n_sig*Del_da[iteff];
      for(size_t itint=0; itint<tint_grid.size(); itint++){
        double t_int = tint_grid[itint];
        if(t_int<t_eff) continue;
        if(t_int > 0.5*T*60.*60.)  continue; // tau_int << T [secs/hours]
        double dX_prev = dX[iT][itint];
        if(da<dX_prev || dX_prev==0) dX[iT][itint] = da;
      }
    }
  }

  //Write out the constraints to text file (function of T and tau_int)
  std::string ofname = "dAlpha"+olabel+".txt";
  std::cout<<"Writing out results (Function of T, tau_int) to: "
    <<ofname<<"\n";
  std::ofstream ofs(ofname);
  ofs<<"tau_int\\T ";
  ofs.precision(3);
  for(auto i_T=0ul; i_T<T_grid.size(); i_T++){
    ofs<<T_grid[i_T]<<" ";
  }
  ofs<<"\n";
  ofs.precision(3);
  for(auto i_tint=0u; i_tint<tint_grid.size(); i_tint++){
    double tau_int = tint_grid[i_tint];
    ofs<<tau_int<<" ";
    for(auto i_T=0u; i_T<T_grid.size(); i_T++){
      ofs<<dX[i_T][i_tint]<<" ";
    }
    ofs<<"\n";
  }
  ofs.close();
}


//******************************************************************************
void outputR_teff(const std::vector<int> &jeff_grid, int tau_avg,
  const std::vector<double> &R_max,
  const std::vector<double> &da_max, const std::vector<double> &Del_da,
  const std::vector<double> &T_obs, const std::string &olabel)
/*
Writes out results (R_max, da_max etc.)
as a function of tau_eff [not tau_int].
Each point will have different T_obs! (also outputted)
*/
{
  // Form:  tau_eff | R da Del_da T_obs
  std::string ofname = "R_teff"+olabel+".txt";
  std::cout<<"Writing out results (Function of tau_eff) to: "
    <<ofname<<"\n";
  std::ofstream ofs(ofname);
  ofs<<"tau_eff R damax Delda Tobs\n";
  ofs.precision(4);

  for(size_t i=0; i<jeff_grid.size(); i++){
    ofs<<jeff_grid[i]*tau_avg<<" "
      <<R_max[i]<<" "
      <<da_max[i]<<" "
      <<Del_da[i]<<" "
      <<T_obs[i]<<"\n";
  }
  ofs.close();

}
