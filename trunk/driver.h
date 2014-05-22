#ifndef D2P_DRIVER_H
#define D2P_DRIVER_H

#include "atom.h"
#include "fix_phonon.h"
#include <list>
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

  void help();
};

#endif
