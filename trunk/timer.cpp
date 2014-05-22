#include "timer.h"

Timer::Timer()
{
  flag = 0;
  start();
return;
}

void Timer::start()
{
 t0 = clock();
 time(&tbeg); tprev = tnow = tbeg;

 flag |= 1;

return;
}

void Timer::stop()
{
  if ( flag&1 ) {
    t2 = clock();
    time(&tnow); tprev = tnow;

    flag |= 2;
  }
return;
}

void Timer::print()
{
  if ( (flag&3) != 3) return;

  cpu_time_used = ((double) (t2 - t0)) / CLOCKS_PER_SEC;
  printf("Total CPU time used: %g seconds; walltime: %g seconds.\n", cpu_time_used, difftime(tnow,tbeg));

return;
}

double Timer::cpu_time()
{
  if ( (flag&3) != 3) return 0.;
  else return ((double) (t2 - t0)) / CLOCKS_PER_SEC;
}

double Timer::wall_time()
{
  if ( (flag&3) != 3) return 0.;
  else return difftime(tnow,tbeg);
}

double Timer::up2now()
{
  if ( (flag&1) != 1) return 0.;
  else {
    time(&tnow);
    return difftime(tnow, tbeg);
  }
}

double Timer::sincelast()
{
  if ( (flag&1) != 1) return 0.;
  else {
    time(&tnow);
    double tinv = difftime(tnow, tprev);
    tprev = tnow;
    return tinv;
  }
}
