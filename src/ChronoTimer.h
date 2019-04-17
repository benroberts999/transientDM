#pragma once
#include <chrono>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <string>
/*
Class that uses std::chrono to easily time code
Usage:
  Create object. By default will start the timer, unless 'true'
  is given as optional input, e.g.:
    ChronoTimer sw(false); //will not start the timing.
    ChronoTimer sw;        //will start

  sw.start() -- starts timing
  sw.stop()  -- 'pauses' timing
  sw.reset() -- clears timer. Needs to be re-stared. Forgets all times.

  start/stop lets you time indevidual sections of code, while not timing
  others.
  Calling start() again will start a new "lap", and save current to total.

  reading_ms() -- returns total elapsed time as double, in ms
  lap_reading_ms() -- same, but only returns time since last start()

  reading_str() and lap_reading_str() -- As above, but outputs as
  formatted string, in units of ms,s,mins, or hours
  (e.g., "1.56 s" or "2.10 hours")

*/

class ChronoTimer {
public:
  ChronoTimer(bool auto_start = true);
  void start();
  void stop();
  void reset();

  double reading_ms() const;
  double lap_reading_ms() const;
  std::string reading_str() const;
  std::string lap_reading_str() const;

private:
  bool running;
  std::chrono::high_resolution_clock::time_point tstart;
  std::string convertHR(double t) const;
  double total_time_ms;
};

//******************************************************************************
inline ChronoTimer::ChronoTimer(bool auto_start)
    : running(false), total_time_ms(0) {
  if (auto_start)
    start(); // true by default
}

//******************************************************************************
inline void ChronoTimer::start() {
  // note: will over-ride any existing reading!
  if (running)
    stop();
  running = true;
  tstart = std::chrono::high_resolution_clock::now();
}

//******************************************************************************
inline void ChronoTimer::stop() {
  if (!running)
    return;

  double current_time = lap_reading_ms();
  running = false; //"turn off" stopwatch

  // update total time
  total_time_ms += current_time;
}

//******************************************************************************
inline void ChronoTimer::reset() {
  running = false;
  total_time_ms = 0;
}

//******************************************************************************
inline double ChronoTimer::lap_reading_ms() const
// Returns value for current riming run (lap)
// Returns double (milliseconds)
{

  if (!running)
    return 0;

  std::chrono::high_resolution_clock::time_point tcurrent =
      std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(tcurrent - tstart)
          .count();

  return ((double)duration) * 1.e-3;
}

//******************************************************************************
inline double ChronoTimer::reading_ms() const
// Returns total value for timnig run
// Returns double (milliseconds)
{
  return lap_reading_ms() + total_time_ms;
}
//******************************************************************************

inline std::string ChronoTimer::reading_str() const {
  return convertHR(reading_ms());
}
//******************************************************************************
inline std::string ChronoTimer::lap_reading_str() const {
  return convertHR(lap_reading_ms());
}

//******************************************************************************
inline std::string ChronoTimer::convertHR(double t) const
// Convers double (in ms) into formmated 2 d.p. string in units of either
// ms, s, mins, or hours, depending on size.
{
  double ot;

  std::string un;
  if (t < 1000) {
    ot = t;
    un = "ms";
  } else if (t < 60000) {
    ot = t / 1000.;
    un = "s";
  } else if (t < 3600000) {
    ot = t / 60000.;
    un = "mins";
  } else {
    ot = t / 3600000.;
    un = "hours";
  }

  std::stringstream ss;
  ss << std::fixed << std::setprecision(2) << ot;
  return ss.str() + " " + un;
}
