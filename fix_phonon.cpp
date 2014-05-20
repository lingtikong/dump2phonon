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

#include "string.h"
#include "fix_phonon.h"
#include "memory.h"

#define MAXLINE 512

/* ---------------------------------------------------------------------- */

FixPhonon::FixPhonon(const char *control)
{
  // initialize allocatables
  nskip = 0;
  memory = NULL;
  fft_data = NULL;
  M_inv_sqrt = NULL;

  Rnow = Rsum = NULL;
  Rqnow = Rqsum = Phi_q = NULL;
  basis = NULL;
  basetype = NULL;
  one = NULL;
  flog = NULL;
  ffw = NULL;

  // initialize default values
  boltz = -1.;
  status = sysdim = 0;
  temperature = 0.;
  type2mass.clear();
  mapfile = prefix = NULL;

  // open control file
  FILE *fp = fopen(control, "r");
  if (fp == NULL){
    printf("\nError: Failed to open file %s!\n", control);
    status = 1;

    return;
  }

  // read the control file
  char str[MAXLINE], units[MAXLINE];
  strcpy(units, "lj");

  fgets(str, MAXLINE, fp);
  while (! feof(fp) ){
    char *ptr = strtok(str, " \n\t\r\f");
    // units
    if (strcmp(ptr, "units") == 0){
      ptr = strtok(NULL, " \n\t\r\f");
      if (strcmp(ptr, "lj") == 0)         boltz = 1.0;
      else if (strcmp(ptr, "real") == 0)  boltz = 0.0019872067;
      else if (strcmp(ptr, "metal") == 0) boltz = 8.617343e-5;
      else if (strcmp(ptr, "si") == 0)    boltz = 1.3806504e-23;
      else if (strcmp(ptr, "cgs") == 0)   boltz = 1.3806504e-16;
      else if (strcmp(ptr, "electron") == 0) boltz = 3.16681534e-6;
      else if (strcmp(ptr, "micro") == 0) boltz = 1.3806504e-8;
      else if (strcmp(ptr, "nano") == 0)  boltz = 0.013806503;
      else {
        printf("\nError: Wrong input, unknown units of: %s\n", ptr);
        status = 2;
        return;
      }
      strcpy(units, ptr);

    // sysdim
    } else if (strcmp(ptr, "sysdim") == 0){
      ptr = strtok(NULL, " \n\t\r\f");
      sysdim = atoi(ptr);
      if (sysdim < 1 || sysdim > 3){
        printf("\nError: Wrong input, sysdim must be within [1, 3].\n");
        status = 2;
        return;
      }

    // temperature
    } else if (strcmp(ptr, "temperature") == 0){
      ptr = strtok(NULL, " \n\t\r\f");
      temperature = atof(ptr);
      if (temperature <= 0.){
        printf("\nError: Wrong input, temperature must be > 0.\n");
        status = 2;
        return;
      }

    // mass
    } else if (strcmp(ptr, "mass") == 0){
      int ip = atoi(strtok(NULL, " \n\t\r\f"));
      double m = atof(strtok(NULL, " \n\t\r\f"));

      if (ip < 1){
        printf("\nError: Wrong input, atomic type must be > 0.\n");
        status = 2;
        return;
      }
      if (m <= 0.){
        printf("\nError: Wrong input, atomic mass must be > 0.\n");
        status = 2;
        return;
      }
      type2mass[ip] = m;

    // map file
    } else if (strcmp(ptr, "map") == 0){
      ptr = strtok(NULL, " \n\t\r\f");
      if (mapfile) delete []mapfile;
      mapfile = new char[strlen(ptr)+1];
      strcpy(mapfile, ptr);

    // prefix for output
    } else if (strcmp(ptr, "prefix") == 0){
      ptr = strtok(NULL, " \n\t\r\f");
      if (prefix) delete []prefix;
      prefix = new char[strlen(ptr)+1];
      strcpy(prefix, ptr);

    // number of frames that should be skipped
    } else if (strcmp(ptr, "nskip") == 0){
      nskip = atoi(strtok(NULL, " \n\t\r\f"));
    }

    fgets(str, MAXLINE, fp);
  }
  fclose(fp);

  // check the completeless of control parameters
  // by default, lj unit
  if (boltz < 0.) boltz = 1.;

  // by default, 3D
  if (sysdim < 1) sysdim = 3;

  // temperature must be supplied
  if (temperature <= 0.){
    printf("\nError: temperature not provided.\n");
    status = 3;
  }
  // by default, map.in
  if (mapfile == NULL){
    mapfile = new char[7];
    strcpy(mapfile, "map.in");
  }
  // by default, phonon
  if (prefix == NULL){
    prefix = new char [7];
    strcpy(prefix, "phonon");
  }

  // open log file, and write some log info
  logfile = new char[strlen(prefix)+5];
  sprintf(logfile,"%s.log", prefix);
  flog = fopen(logfile, "w");
  fprintf(flog, "# Computation of dynamical matrix from MD trajectories.\n");
  fprintf(flog, "# Control parameters read from: %s\n", control);
  fprintf(flog, "# Key parameter info:\n");
  fprintf(flog, "    units        %s\n", units);
  fprintf(flog, "    sysdim       %d\n", sysdim);
  fprintf(flog, "    temperature  %g\n", temperature);
  fprintf(flog, "    mapfile      %s\n", mapfile);
  fprintf(flog, "    prefix       %s\n", prefix);
  for (std::map<int, double>::iterator it = type2mass.begin(); it != type2mass.end(); ++it){
    fprintf(flog, "    mass     %d   %g\n", it->first, it->second);
  }
  fprintf(flog, "# Boltzmann constant under current units: %lg\n", boltz);
  fflush(flog);
return; 
}

/* ---------------------------------------------------------------------- * 
 * Private method to initialize FixPhonon
 * Before calling this method, one frame of the trajectory must be read and
 * passed to "one"
 * ---------------------------------------------------------------------- */
int FixPhonon::init()
{
  memory = new Memory();

  // mapping index
  tag2surf.clear(); // clear map info
  surf2tag.clear();

  // get the mapping between lattice indices and atom IDs
  if (readmap()) return 1;

  // create FFT and allocate memory for FFT
  npt = nq = nx*ny*nz;

  fft_dim   = nucell*sysdim;
  fft_dim2  = fft_dim*fft_dim;

  memory->create(fft_data, nq, "fix_phonon:fft_data");
  ffw = fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(fft_data), reinterpret_cast<fftw_complex*>(fft_data), FFTW_FORWARD, FFTW_MEASURE);

  // allocate variables
  memory->create(Rnow, npt, fft_dim, "fix_phonon:Rnow");
  memory->create(Rsum, npt, fft_dim, "fix_phonon:Rsum");
                              
  memory->create(basis, nucell, sysdim, "fix_phonon:basis");

  // because of hermit, only nearly half of q points are stored
  memory->create(Rqnow, nq, fft_dim,  "fix_phonon:Rqnow");
  memory->create(Rqsum, nq, fft_dim2, "fix_phonon:Rqsum");
  memory->create(Phi_q, nq, fft_dim2, "fix_phonon:Phi_q");

  // output some information on the system to log file
  for (int i = 0; i < 15; ++i) fprintf(flog,"####"); fprintf(flog,"\n");
  fprintf(flog,"# total number of atoms in the group       : %d\n", ngroup);
  fprintf(flog,"# dimension of the system                  : %d D\n", sysdim);
  fprintf(flog,"# number of atoms per unit cell            : %d\n", nucell);
  fprintf(flog,"# dimension of the FFT mesh                : %d x %d x %d\n", nx, ny, nz);
  for (int i = 0; i < 60; ++i) fprintf(flog,"#"); fprintf(flog,"\n");
  fprintf(flog,"# mapping information between lattice indx and atom id\n");
  fprintf(flog,"# nx ny nz nucell\n");
  fprintf(flog,"%d %d %d %d\n", nx, ny, nz, nucell);
  fprintf(flog,"# l1 l2 l3 k atom_id\n");
  int ix, iy, iz, iu;
  for (int idx = 0; idx < ngroup; ++idx){
    int id = surf2tag[idx];
    iu   = idx%nucell;
    iz   = (idx/nucell)%nz;
    iy   = (idx/(nucell*nz))%ny;
    ix   = (idx/(nucell*nz*ny))%nx;
    fprintf(flog,"%d %d %d %d %d\n", ix, iy, iz, iu, id);
  }
  for (int i = 0; i < 15; ++i) fprintf(flog,"####"); fprintf(flog,"\n");
  fflush(flog);
 
  // initialize mass and type info
  memory->create(M_inv_sqrt, nucell, "M_inv_sqrt");
  memory->create(basetype, nucell, "basetype");
  for (int iu = 0; iu < nucell; ++iu){
    int id = surf2tag[iu];
    int ip = one->attyp[id];
    if (type2mass.count(ip) < 1){
      printf("\nError: mass info for type %d not defined!\n", ip);
      return 2;
    }
    M_inv_sqrt[iu] = 1./sqrt(type2mass[ip]);
    basetype[iu]   = ip;
  }

return 0;
} // end of init

FixPhonon::~FixPhonon()
{
  // delete locally stored array
  memory->destroy(Rnow);
  memory->destroy(Rsum);

  memory->destroy(basis);

  memory->destroy(Rqnow);
  memory->destroy(Rqsum);
  memory->destroy(Phi_q);

  if (prefix) delete []prefix;
  if (mapfile) delete []mapfile;
  if (logfile) delete []logfile;

  memory->destroy(M_inv_sqrt);
  memory->destroy(basetype);

  // destroy FFT
  memory->destroy(fft_data);
  if (ffw) fftw_destroy_plan(ffw);
  
  // clear map info
  tag2surf.clear();
  surf2tag.clear();

  if (flog) fclose(flog);

return;
}

/* ---------------------------------------------------------------------- */

void FixPhonon::setup()
{
  // initialize accumulating variables
  for (int i = 0; i < npt; ++i)
  for (int j = 0; j < fft_dim;  ++j) Rsum[i][j] = 0.;

  for (int i = 0; i < nq; ++i)
  for (int j = 0; j < fft_dim2; ++j) Rqsum[i][j] = std::complex<double> (0.,0.);

  for (int i = 0; i < 6; ++i) hsum[i] = 0.;

  for (int i = 0; i < nucell; ++i)
  for (int j = 0; j < sysdim; ++j) basis[i][j] = 0.;

  neval  = 0;
return;
}

/* ---------------------------------------------------------------------- */

void FixPhonon::end_of_step()
{
  // get Rnow
  for (int il = 0; il < npt; ++il)
  for (int iu = 0; iu < nucell; ++iu){
    int idx = il * nucell + iu;
    int id = surf2tag[idx];
    int ndim = iu * sysdim;
    for (int idim = 0; idim < sysdim; ++idim) Rnow[il][ndim+idim] = one->atpos[id][idim];
  }

  // get Rsum
  for (int il = 0; il < npt;  ++il)
  for (int idim = 0; idim < fft_dim; ++idim) Rsum[il][idim] += Rnow[il][idim];

  // FFT R(r) to get R(q)
  for (int idim = 0; idim < fft_dim; ++idim){
    for (int il = 0; il < npt; ++il) fft_data[il] = std::complex<double>(Rnow[il][idim], 0.);

    fftw_execute(ffw);

    for (int idq = 0; idq < nq; ++idq) Rqnow[idq][idim] = fft_data[idq];
  }

  // to get sum(R(q).R(q)*)
  for (int idq = 0; idq < nq; ++idq){
    int ndim = 0;
    for (int idim = 0; idim < fft_dim; ++idim)
    for (int jdim = 0; jdim < fft_dim; ++jdim) Rqsum[idq][ndim++] += Rqnow[idq][idim]*conj(Rqnow[idq][jdim]);
  }

  // get basis info
  if (fft_dim > sysdim){
    double d2o[3];
    for (int il = 0; il < npt; ++il){
      for (int iu = 1; iu < nucell; ++iu){
        int i0 = surf2tag[il*nucell];
        int id = surf2tag[il*nucell+iu];
        for (int idim = 0; idim < sysdim; ++idim) d2o[idim] = one->atpos[id][idim] - one->atpos[i0][idim];
        one->apply_pbc(d2o);
        for (int idim = 0; idim < sysdim; ++idim) basis[iu][idim] += d2o[idim];
      }
    }
  }

  // get lattice vector info
  for (int i = 0; i < 6; ++i) hsum[i] += one->h[i];

  // increment counter
  ++neval;

return;
}   // end of end_of_step()

/* ----------------------------------------------------------------------
 * private method, to read the mapping info from file
 * --------------------------------------------------------------------*/
int FixPhonon::readmap()
{
  // auto-generate mapfile for "cluster" (gamma only system)
  if (strcmp(mapfile, "GAMMA") == 0){
    nx = ny = nz = ntotal = 1;
    ngroup = nucell = one->natom;

    for (int id = 1; id <= one->natom; ++id){
      tag2surf[id] = id-1;
      surf2tag[id-1] = id;
    }

    return 0;
  }

  // read from map file
  char line[MAXLINE];
  FILE *fp = fopen(mapfile, "r");
  if (fp == NULL){
    printf("\nError: Cannot open input map file: %s!\n", mapfile);

    return 1;
  }

  if (fgets(line,MAXLINE,fp) == NULL){
    printf("\nError while reading header of the map file!\n");

    return 2;
  }
  nx     = atoi(strtok(line, " \n\t\r\f"));
  ny     = atoi(strtok(NULL, " \n\t\r\f"));
  nz     = atoi(strtok(NULL, " \n\t\r\f"));
  nucell = atoi(strtok(NULL, " \n\t\r\f"));
  ntotal = nx*ny*nz;
  ngroup = ntotal*nucell;
  if (ngroup > one->natom){
    printf("\nError: %d atoms found in the trajectory file, while %d suggested by the map file!\n", ngroup, one->natom);

    return 3;
  }
  
  // second line of mapfile is comment
  if (fgets(line,MAXLINE,fp) == NULL){
    printf("\nError while reading comment of the map file!");

    return 3;
  }

  int info = 0;
  // the remaining lines carry the mapping info
  for (int i = 0; i < ngroup; ++i){
    if (fgets(line,MAXLINE,fp) == NULL) {info = 1; break;} 
    int ix = atoi(strtok(line, " \n\t\r\f"));
    int iy = atoi(strtok(NULL, " \n\t\r\f"));
    int iz = atoi(strtok(NULL, " \n\t\r\f"));
    int iu = atoi(strtok(NULL, " \n\t\r\f"));
    int id = atoi(strtok(NULL, " \n\t\r\f"));

    // check if index is in correct range
    if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz|| iu < 0 || iu >= nucell) {info = 2; break;}
    // 1 <= id <= natoms
    if (id < 1 || id > one->natom) {info = 3; break;}
    int idx = ((ix*ny+iy)*nz+iz)*nucell + iu;
    tag2surf[id]  = idx;
    surf2tag[idx] = id;
  }
  fclose(fp);

  if (info > 0){
    printf("\nError while reading the map file %s!\n", mapfile);
    return 4;
  }
  if (tag2surf.size() != surf2tag.size() || tag2surf.size() != ngroup){
    printf("\nError: The mapping info read from %s is incomplete!\n", mapfile);

    return 4;
  }
  
return 0;
}

/* ----------------------------------------------------------------------
 * private method, to output the force constant matrix
 * --------------------------------------------------------------------*/
void FixPhonon::post_run()
{
  if (neval < 1) return;
  fprintf(flog, "# %d frames were used.\n", neval);

  double inv_neval = 1. /double(neval);

  // to get <Rq.Rq*>
  for (int idq = 0; idq < nq; ++idq)
  for (int idim = 0; idim < fft_dim2; ++idim) Phi_q[idq][idim] = Rqsum[idq][idim] * inv_neval;

  // to get <R>
  for (int il = 0; il < npt; ++il)
  for (int idim = 0; idim < fft_dim; ++idim) Rnow[il][idim] = Rsum[il][idim] * inv_neval;

  // to get <R>q
  for (int idim = 0; idim < fft_dim; ++idim){
    for (int il = 0; il < npt; ++il) fft_data[il] = std::complex<double>(Rnow[il][idim], 0.);

    fftw_execute(ffw);

    for (int idq = 0; idq < nq; ++idq) Rqnow[idq][idim]  = fft_data[idq];
  }

  // to get G(q) = <Rq.Rq*> - <R>q.<R*>q
  for (int idq = 0; idq < nq; ++idq){
    int ndim = 0;
    for (int idim = 0; idim < fft_dim; ++idim)
    for (int jdim = 0; jdim < fft_dim; ++jdim) Phi_q[idq][ndim++] -= Rqnow[idq][idim]*conj(Rqnow[idq][jdim]);
  }

  // to get Phi = KT.G^-1; normalization of FFTW data is done here
  double kBT = boltz * temperature * double(nq);
  for (int idq = 0; idq < nq; ++idq){
    GaussJordan(fft_dim, Phi_q[idq]);

    for (int idim = 0; idim < fft_dim2; ++idim) Phi_q[idq][idim] *= kBT;
  }

  // get basis info
  for (int iu = 1; iu < nucell; ++iu)
  for (int idim = 0; idim < sysdim; ++idim) basis[iu][idim] /= double(npt)*double(neval);

  // get unit cell base vector info; might be incorrect if MD pbc and FixPhonon pbc mismatch.
  double basevec[9];
  basevec[1] = basevec[2] = basevec[5] = 0.;
  basevec[0] = hsum[0] * inv_neval / double(nx);
  basevec[4] = hsum[1] * inv_neval / double(ny);
  basevec[8] = hsum[2] * inv_neval / double(nz);
  basevec[7] = hsum[3] * inv_neval / double(nz);
  basevec[6] = hsum[4] * inv_neval / double(nz);
  basevec[3] = hsum[5] * inv_neval / double(ny);
    
  // write binary file, in fact, it is the force constants matrix that is written
  // Enforcement of ASR and the conversion of dynamical matrix is done in the postprocessing code
  char fname[MAXLINE];
  sprintf(fname,"%s.bin", prefix);
  FILE *fp_bin = fopen(fname,"wb");

  fwrite(&sysdim, sizeof(int),    1, fp_bin);
  fwrite(&nx,     sizeof(int),    1, fp_bin);
  fwrite(&ny,     sizeof(int),    1, fp_bin);
  fwrite(&nz,     sizeof(int),    1, fp_bin);
  fwrite(&nucell, sizeof(int),    1, fp_bin);
  fwrite(&boltz,  sizeof(double), 1, fp_bin);

  fwrite(Phi_q[0],sizeof(double), nq*fft_dim2*2, fp_bin);

  fwrite(&temperature,  sizeof(double),1,      fp_bin);
  fwrite(&basevec[0],   sizeof(double),9,      fp_bin);
  fwrite(basis[0],      sizeof(double),fft_dim,fp_bin);
  fwrite(basetype,      sizeof(int),   nucell, fp_bin);
  fwrite(M_inv_sqrt,    sizeof(double),nucell, fp_bin);

  fclose(fp_bin);

return;
}   // end of postprocess

/* ----------------------------------------------------------------------
 * private method, to get the inverse of a complex matrix by means of
 * Gaussian-Jordan Elimination with full pivoting; square matrix required.
 *
 * Adapted from the Numerical Recipes in Fortran.
 * --------------------------------------------------------------------*/
void FixPhonon::GaussJordan(int n, std::complex<double> *Mat)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
  int *indxc,*indxr,*ipiv;
  double big, nmjk;
  std::complex<double> dum, pivinv;

  indxc = new int[n];
  indxr = new int[n];
  ipiv  = new int[n];

  for (i = 0; i < n; ++i) ipiv[i] = 0;
  for (i = 0; i < n; ++i){
    big = 0.;
    for (j = 0; j < n; ++j){
      if (ipiv[j] != 1){
        for (k = 0; k < n; ++k){
          if (ipiv[k] == 0){
            idr = j*n+k;
            nmjk = norm(Mat[idr]);
            if (nmjk >= big){
              big  = nmjk;
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1){
            printf("\nError: Singular matrix in complex GaussJordan!");
            exit(1);
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l = 0; l < n; ++l){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == std::complex<double>(0.,0.)){
      printf("\nError: Singular matrix in complex GaussJordan!");
      exit(1);
    }
    
    pivinv = 1./ Mat[idr];
    Mat[idr] = std::complex<double>(1.,0.);
    idr = icol*n;
    for (l = 0; l < n; ++l) Mat[idr+l] *= pivinv;
    for (ll = 0; ll < n; ++ll){
      if (ll != icol){
        idc = ll*n + icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l = 0; l < n; ++l) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }

  for (l = n-1; l >= 0; --l){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k = 0; k < n; ++k){
        idr = k*n + rl;
        idc = k*n + cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }
  delete []indxr;
  delete []indxc;
  delete []ipiv;

return;
}

/* --------------------------------------------------------------------*/
