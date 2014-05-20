#ifndef D2P_DRIVER_H
#define D2P_DRIVER_H

#include "atom.h"
#include "memory.h"
#include "fix_phonon.h"
#include <list>
//#include <vector>
#include <string>

class Driver {
public:
  Driver(int, char**);
  ~Driver();

private:
  std::string fname;
  std::list<std::string> flist;

  DumpAtom *one, *ref;                    // pointer to one frame

  int first, iframe;
  FILE *fdump;
  void readdump();
  void compute();

  void help();
};

#endif
