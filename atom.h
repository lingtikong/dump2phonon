#ifndef DUMP_ATOM_H
#define DUMP_ATOM_H

/* -----------------------------------------------------------------------------
 * A class which contains one frame of the atom style dump of lammps
 * -------------------------------------------------------------------------- */

#include "memory.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <map>
#include <string>
#include "global.h"

using namespace std;

class DumpAtom {
public:
  DumpAtom(FILE *fp, const char *);
  ~DumpAtom();

  int iframe;
  int natom, ntype, tstep;
  int initialized;

  char *fname;

  Memory *memory;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz;
  double lx, ly, lz, h_inv[6], h[6];
  double hx, hy, hz;
  int *attyp;     // note: atom IDs go from 1 to natom; type ID from 1 to ntype
  double **atpos;
  double axis[3][3];

  void wrap2ref(DumpAtom *);
  void apply_pbc(double *);

private:
  int cartesian, triclinic;
  double **x, **s;
  void car2dir();
  void dir2car();

  int count_words(const char *);
};
#endif
