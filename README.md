# Fibre-network transient DM search

Compile using Makefile: (run _$ make_ from command line).
The relevant executable is called _transientDM_

Program reads input parameters from a text file: _transientDM.in_.
Another text file (by default called: _clocklist.in_) contains a list of
data files that will be read in by the program.

Another program, _combinedX_ (reads input from _combinedX.in_) takes in a bunch of given data files (outputs of _transientDM_), and forms the resulting constraints [on delta(alpha)(tau_int), Lambda(T), and Lambda(d)]


## Input options description
Each relevent option is explained inside the _transientDM.in_. Here is a basic list of what you can specify:
  * For which time period to run the analysis
  * The 'averaging time' (in seconds) -- default is 60s, will average data into 60s chunks before the analysis. By default, it will throw away any data that is not perfectly continuous over this averaging period, but you can specify how many missing/skipped points are acceptable (0 by default). *See note (1)*
  * Which tau_int's (interaction/transient duration) to look for explicitly. Minimum must be >= averaging time.
  * Minimum number of clocks to use, and whether or not to 'force' the code to use certain clock combos. If these are set to true, program will ignore any time periods that these specific clocks were not active. By choosing only the best clocks, we can get a much improved constraint (but it applies to only a much smaller observation time)
  * Use Gaussian or flat profile. Gaussian is default. Makes essentially no difference.
  * n_sig sets the confidence level, limit = |max| + n_sig * error
  (n=1 => 68% CL, n=1.645 => 90% etc.)


### notes
 (1) Note: method assumes white noise, so the data must first be averaged
 to larger than the worst servo time (~60s).


## some example code:

```cpp
int main() {
  int x = 5 + 6;
  ClockNetwork net(filenames,tau_avg,max_bad_avg);
  net.replaceWithRandomNoise();
}
```
