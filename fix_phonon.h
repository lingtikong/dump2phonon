/* ----------------------------------------------------------------------
   Contributing authors:
     Ling-Ti Kong

   Contact:
     School of Materials Science and Engineering,
     Shanghai Jiao Tong University,
     800 Dongchuan Road, Minhang,
     Shanghai 200240, CHINA

     konglt@sjtu.edu.cn; konglt@gmail.com
------------------------------------------------------------------------- */
#ifndef FIXPHONON_H
#define FIXPHONON_H

#include <map>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include <complex>
#include "fftw3.h"
#include "memory.h"

class FixPhonon {
 public:
  FixPhonon(const char *);
  ~FixPhonon();

  int init();
  void setup();
  void end_of_step();
  void post_run();

  int status;                                   // 0, no problem; > 0 error at various levels.
  int nskip;
  DumpAtom *one;

 private:
  Memory *memory;

  int nx,ny,nz,nucell,ntotal;                   // surface dimensions in x- and y-direction, number of atom per unit surface cell
  int neval;                                    // # of evaluations
  int sysdim;                                   // system dimension
  int ngroup;                                   // total number of atoms in group
  FILE *flog;
  char *prefix, *logfile;
  
  std::map<int,double> type2mass;
  double *M_inv_sqrt;
  double boltz;

  fftw_plan  ffw;                               // to do fft via the fft3d wraper
  std::complex<double> *fft_data;
  int npt, nq;
  int fft_dim, fft_dim2;
  
  std::map<int,int> tag2surf, surf2tag;         // Mapping info

  double **Rnow;                                // Current R(r) on local proc
  double **Rsum;                                // Accumulated R(r) on local proc

  int *recvcnts, *displs;                       // MPI related variables

  std::complex<double> **Rqnow;                 // Current R(q) on local proc
  std::complex<double> **Rqsum;                 // Accumulator for conj(R(q)_alpha)*R(q)_beta
  std::complex<double> **Phi_q;                 // Phi's on local proc

  int readmap();                                // to read the mapping of gf atoms
  char *mapfile;                                // file name of the map file

  double temperature;
  double hsum[6], **basis;
  int *basetype;

  // private methods to do matrix inversion
  void GaussJordan(int, std::complex<double>*);

};
#endif
