# Fibre-network transient DM search

Compile using Makefile: (run _$ make_ from command line).
The relevant executable is called _transientDM_

Program reads input parameters from a text file: _transientDM.in_.
Another text file (by default called: _clocklist.in_) contains a list of
data files that will be read in by the program.
  * data files are assumed to be in format 'PREFIX_CLOCK1-CLOCK2.SUFFIX'
  * prefix + suffix are input options

Another program, _combinedX_ (reads input from _combinedX.in_) takes in a bunch of given data files (outputs of _transientDM_), and forms the resulting constraints [on delta(alpha)(tau_int), Lambda(T), and Lambda(d)]

**Note:** you must first create a directory called _/obj/_ (the compiled, but un-linked) programs go here - the reason is that then only parts of the programs that actually change need to be re-compiles.
Technically, these files are not needed after compilation, so feel free to delete them.


## Input options description
Each relevant option is explained inside the _transientDM.in_. Here is a basic list of what you can specify:
  * For which time period to run the analysis
  * The 'averaging time' (in seconds) -- default is 60s, will average data into 60s chunks before the analysis. By default, it will throw away any data that is not perfectly continuous over this averaging period, but you can specify how many missing/skipped points are acceptable (0 by default). *See note (1)*
  * Which tau_int's (interaction/transient duration) to look for explicitly. Minimum must be >= averaging time.
  * Minimum number of clocks to use, and whether or not to 'force' the code to use certain clock combos. If these are set to true, program will ignore any time periods that these specific clocks were not active. By choosing only the best clocks, we can get a much improved constraint (but it applies to only a much smaller observation time)
  * Use Gaussian or flat profile. Gaussian is default. Makes essentially no difference.
  * n_sig sets the confidence level, limit = |max| + n_sig * error
  (n=1 => 68% CL, n=1.645 => 90% etc.)
  * Option to inject a fake event (either into real data, or randomly generated noise). This gives a good test of the method.


### notes
 (1) Note: method assumes white noise, so the data must first be averaged
 to larger than the worst servo time (~60s).


## Files (all located in _/src/_):

  * **transientDM**  -  Main program. Contains main() and a few other data analysis routines (e.g., calculate the constraints, and the threshold etc.).
  * **combinedX**  -  Combines several output text files (output from transientDM) that contain constraints as a function of tau_int and T, and combines them into single constrains file (in same form)
  * **ClockInfo.h**  -  Header that contains details for each clock (location, sensitivity coeficients etc.) Needs to be updated if new clocks are introduced to the network
  * **ClockNetwork**  -  Class that holds the clock data, and does all the actual data processing (calculating likelihoods, standard deviations, averages the data etc.). It also contains the routine that parses the data files into the relevant arrays. Minor changes might be needed here if the file format of the data files changes. Note: For now, explicitly assumes the clock network is essentially co-located so the signals are coincident. It's fairly simple to extend this to allow for more general case.
  * **DMs_signalTemplates**  -  Contains functions that define the transient DM signal in the clock data. Has options for either Gaussian-profile DM objects, or flat-profile (top-hat profile).
  * **DataIO**  -  Functions to read in data from a text file into arrays. At the moment, only reads a two-collumn file. Trivial to extend to 3-col
  * **RNG_randomNumberGenerators**  -  provides some simple thread-safe functions to (relatively) quickly generate random real numbers. Can either generate according to a flat distribution (between given min,max), or with a Gaussian distribution [using inverse transform sampling] (with given mean and standard deviation). Can easily update for other distributions too.
  * **ChronoTimer**  -  Just simple class to help accurately time code. Relies on actual time, not CPU clicks, so can be used to accurately time parallel regions of code.


### some example code:

  *  This shows roughly the fitting procedure for a single value of tau_int and a single interaction time (t0).

```cpp
int main() {

  // Opens the file "clocklist.in" (which contains paths to data files)
  // Reads the data file names/locations into list: filenames
  std::vector<std::string> filenames;
  DataIO::getFileNames(filenames, "/path/to/clocklist.in");

  // Generate "ClockNetwork" object:
  int tau_avg = 60;    // average the data into 60s chunks
  int max_bad_avg = 0; // don't allow any missing data points
  ClockNetwork net(filenames, tau_avg, max_bad_avg);

  // Generate the signal template (template is signal, with deltaAlpha=1)
  std::vector<std::vector<double>> s; // array to hold signal
  int j_int = tau_int / tau_avg;      // Interaction duration (in epochs)
  int nJW = 1;                        // Take 'window' to be 1*tau_int
  net.genSignalTemplate(s, nJW, j_int);

  // Time to check for a DM interaction:
  int j0 = t0 / tau_avg;

  // Form the sub-network of independent clock-pairs:
  std::vector<int> indep_pairs;
  net.formIndependentSubnet(indep_pairs, j0, Jw);

  // Calculate the 'dHs' and 'sHs' values (see paper for definition)
  auto xHs = net.calculate_dHs_sHs(indep_pairs, s, j0);

  // Now, we can get the best-fit \delta\alpha, R, and \Delta\alpha:
  double da_bf = fabs(xHs.dHs) / xHs.sHs;       // best-fit amplitude
  double Delta = 1. / sqrt(xHs.sHs);            // Error term
  double R = fabs(xHs.dHs) / sqrt(2 * xHs.sHs); // SNR ratio

  // Normally, will be looped over all relevant values of tau_int.
  // For each tau_int, loops over all values of t0 (times), maximises R and da
}
```


## Note:

This code uses two assumptions:
  1) All clocks affected by TD simultaneously (good for tau >~ 15s)
  2) All noise in white Gaussian frequency noise

To correct (1), need to update [in ClockNetwork.cpp]:
  * Need to store positions of all clocks (this is done, but not used)
  * genSignalTemplate() routine needs to be updated to take into account
  the position of each clock, and the incident speed/direction of TD
  (see my PRD 2016 paper for formulas)

To correct (2), need to update:
 * Need to calculate (or read in) covariance matrix + inverse
 * Need to update calculate_dHs_sHs() function to use matrix


Also: it is hard-coded to use the Sr/Yb+/Hg clocks of PTB/SYRTE/NPL.
To use other clocks, must add their info to 'ClockInfo.h' header
