#ifndef TIMER_H
#define TIMER_H

#include "stdio.h"
#include "stdlib.h"
#include "time.h"

class Timer {
public:
  Timer();

  void start();       // start the timer
  void stop();        // stop the timer
  void print();       // print the time info

  double cpu_time();  // return cpu time used, if timer stopped
  double wall_time(); // return wall time used, if timer stoped

  double up2now();    // return wall time since start
  double sincelast(); // return wall time since last call of timer

private:
  clock_t t0, t2;
  double cpu_time_used;

  time_t tbeg, tnow, tprev;

  int flag;
};

#endif
