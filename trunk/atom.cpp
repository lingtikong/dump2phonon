#include "atom.h"

/*------------------------------------------------------------------------------
 * Constructor, to read one image from the atom style dump file of lammps
 *------------------------------------------------------------------------------
 * Parameters:
 *  fp       : (in) File pointer, of dump file
 *  dumpfile : (in) dump file
 *----------------------------------------------------------------------------*/
DumpAtom::DumpAtom(FILE *fp, const char *dumpfile)
{
  natom = ntype = tstep = 0;
  initialized = triclinic = 0;
  xy = xz = yz = 0.;
  cartesian = 0;

  attyp = NULL; atpos = x = s = NULL;
  memory = NULL; fname = NULL;

  if (fp == NULL) return;

  fname = new char [strlen(dumpfile)+1];
  strcpy(fname, dumpfile);

  // time step info
  char str[MAXLINE];
  if (fgets(str,MAXLINE, fp) == NULL) return;
  fgets(str,MAXLINE, fp);
  tstep = atoi(strtok(str, " \n\t\r\f"));

  // # of atoms info
  fgets(str,MAXLINE, fp);
  fgets(str,MAXLINE, fp);
  natom = atoi(strtok(str, " \n\t\r\f"));
  if (natom < 1) return;

  // box bounds info
  fgets(str,MAXLINE, fp);
  fgets(str,MAXLINE, fp);
  int n = count_words(str);
  xlo = atof(strtok(str, " \n\t\r\f"));
  xhi = atof(strtok(NULL," \n\t\r\f"));
  if (n == 3) xy = atof(strtok(NULL," \n\t\r\f"));

  fgets(str,MAXLINE, fp);
  n = count_words(str);
  ylo = atof(strtok(str, " \n\t\r\f"));
  yhi = atof(strtok(NULL," \n\t\r\f"));
  if (n == 3) xz = atof(strtok(NULL," \n\t\r\f"));

  fgets(str,MAXLINE, fp);
  n = count_words(str);
  zlo = atof(strtok(str, " \n\t\r\f"));
  zhi = atof(strtok(NULL," \n\t\r\f"));
  if (n == 3) yz = atof(strtok(NULL," \n\t\r\f"));
  
  if (xy*xy+xz*xz+yz*yz > ZERO){
    triclinic = 1;

    xlo -= MIN(MIN(0., xy), MIN(xz, xy+xz));
    xhi -= MAX(MAX(0., xy), MAX(xz, xy+xz));
    ylo -= MIN(0., yz);
    yhi -= MAX(0., yz);
  }

  // fields info
  int dcols[6], fcord = 7;
  for (int i = 0; i < 6; ++i) dcols[i] = i;
  fgets(str,MAXLINE, fp);
  char *ptr = strtok(str, " \n\t\r\f");
  for (int i = 0; i < 2; ++i) ptr = strtok(NULL," \n\t\r\f");
  int ic = 1;
  while (ptr){
    if (strcmp(ptr, "id") == 0)   dcols[1] = ic;
    if (strcmp(ptr, "type") == 0) dcols[2] = ic;
    if (strcmp(ptr, "xs") == 0){  dcols[3] = ic; fcord |= 1;}
    if (strcmp(ptr, "ys") == 0){  dcols[4] = ic; fcord |= 2;}
    if (strcmp(ptr, "zs") == 0){  dcols[5] = ic; fcord |= 4;}
    if (strcmp(ptr, "x") == 0){   dcols[3] = ic; fcord &= 6;}
    if (strcmp(ptr, "y") == 0){   dcols[4] = ic; fcord &= 5;}
    if (strcmp(ptr, "z") == 0){   dcols[5] = ic; fcord &= 3;}

    ++ic;
    ptr = strtok(NULL," \n\t\r\f");
  }
  fcord &= 7;
  if (fcord != 7 && fcord != 0) return;

  memory = new Memory();
  memory->create(attyp, natom+1, "attyp");

  // always assumes fractional coordinate
  memory->create(s, natom+1, 3, "s");
  atpos = s;

  // read coordinate
  int id, ip;
  double xp, yp, zp;
  for (int i = 0; i < natom; ++i){
    fgets(str,MAXLINE, fp);

    int ic = 1, frd = 0;
    ptr = strtok(str, " \n\t\r\f");
    while (ptr){
      if (ic == dcols[1]){ id = atoi(ptr); frd |=  1; }
      if (ic == dcols[2]){ ip = atoi(ptr); frd |=  2; }
      if (ic == dcols[3]){ xp = atof(ptr); frd |=  4; }
      if (ic == dcols[4]){ yp = atof(ptr); frd |=  8; }
      if (ic == dcols[5]){ zp = atof(ptr); frd |= 16; }

      ptr = strtok(NULL," \n\t\r\f"); ++ic;
    }
    if (frd == 31){
      attyp[id] = ip; ntype = MAX(ip, ntype);
      atpos[id][0] = xp;
      atpos[id][1] = yp;
      atpos[id][2] = zp;
    } else { return; } // insufficient info, return
  }

  lx = h[0] = axis[0][0] = xhi - xlo;
  ly = h[1] = axis[1][1] = yhi - ylo;
  lz = h[2] = axis[2][2] = zhi - zlo;
  hx = 0.5*lx;  hy = 0.5*ly;  hz = 0.5*lz;

  h[5] = axis[1][0] = xy;
  h[4] = axis[2][0] = xz;
  h[3] = axis[2][1] = yz;
  axis[0][1] = axis[0][2] = axis[1][2] = 0.;

  h_inv[0] = 1./h[0];
  h_inv[1] = 1./h[1];
  h_inv[2] = 1./h[2];
  h_inv[3] = -h[5] / (h[1]*h[2]);
  h_inv[4] = (h[5]*h[3] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
  h_inv[5] = -h[3] / (h[0]*h[1]);

  // in case cartesian coordinate read, convert to fractional
  if ((fcord&7) == 0){
    x = s;
    s = NULL;
    car2dir();
  }

  initialized = 1;

return;
}

/*------------------------------------------------------------------------------
 * Deconstructor, to free allocated memory
 *----------------------------------------------------------------------------*/
DumpAtom::~DumpAtom()
{
  if (fname) delete []fname;
  if (attyp) memory->destroy(attyp);

  atpos = NULL;
  if (x) memory->destroy(x);
  if (s) memory->destroy(s);

  if (memory) delete memory;
}


/*------------------------------------------------------------------------------
 * Method to convert fractional coordinates into cartesian
 * s must be available
 *----------------------------------------------------------------------------*/
void DumpAtom::dir2car()
{
  if (cartesian) return;

  if (x == NULL){
    memory->create(x, natom+1, 3,"x");

    if (triclinic){
      for (int id = 1; id <= natom; ++id){
        x[id][0] = s[id][0]*lx + s[id][1]*xy + s[id][2]*xz + xlo;
        x[id][1] = s[id][1]*ly + s[id][2]*yz + ylo;
        x[id][2] = s[id][2]*lz + zlo;
      }

    } else {
      for (int id = 1; id <= natom; ++id){
        x[id][0] = s[id][0]*lx + xlo;
        x[id][1] = s[id][1]*ly + ylo;
        x[id][2] = s[id][2]*lz + zlo;
      }
    }
  }

  atpos = x;
  cartesian = 1;
return;
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinates into fractional
 * x must be available
 *----------------------------------------------------------------------------*/
void DumpAtom::car2dir()
{
  if (cartesian == 0) return;

  if (s == NULL){
    memory->create(s, natom+1, 3,"s");
    
    double x0[3];
    if (triclinic){
      for (int id = 1; id <= natom; ++id){
        x0[0] = x[id][0] - xlo;
        x0[1] = x[id][1] - ylo;
        x0[2] = x[id][2] - zlo;
        s[id][0] = h_inv[0]*x0[0] + h_inv[5]*x0[1] + h_inv[4]*x0[2];
        s[id][1] = h_inv[1]*x0[1] + h_inv[3]*x0[2];
        s[id][2] = h_inv[2]*x0[2];
      }

    } else {

      for (int id = 1; id <= natom; ++id){
        s[id][0] = h_inv[0]*(x[id][0] - xlo);
        s[id][1] = h_inv[1]*(x[id][1] - ylo);
        s[id][2] = h_inv[2]*(x[id][2] - zlo);
      }
    }
  }

  atpos = s;
  cartesian = 0;
return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int DumpAtom::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy, n, "copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) ++n;

  memory->sfree(copy);
  return n;
}

/* -----------------------------------------------------------------------------
 * Private method to wrap atoms back to its original/reference image
 * ---------------------------------------------------------------------------*/
void DumpAtom::wrap2ref(DumpAtom *ref)
{
  if (natom != ref->natom){
    printf("\nError: # of atoms from %s-%d and %s-%d mismatch!\n", fname, iframe, ref->fname, ref->iframe);

    return;
  }

  // need cartesian coordinate system
  dir2car();
  ref->dir2car();

  if (triclinic) {  // non-orthogonal box
    for (int id = 1; id <= natom; ++id){
      double xij = atpos[id][0] - ref->atpos[id][0];
      double yij = atpos[id][1] - ref->atpos[id][1];
      double zij = atpos[id][2] - ref->atpos[id][2];

      while (zij >= hz){
        xij -= xz;
        yij -= yz;
        zij -= lz;
      }
      while (zij < -hz){
        xij += xz;
        yij += yz;
        zij += lz;
      }
      
      while (yij >= hy){
        xij -= xy;
        yij -= ly;
      }
      while (yij < -hy){
        xij += xy;
        yij += ly;
      }
      
      while (xij >= hx) xij -= lx;
      while (xij < -hx) xij += lx;
   
      atpos[id][0] = ref->atpos[id][0] + xij;
      atpos[id][1] = ref->atpos[id][1] + yij;
      atpos[id][2] = ref->atpos[id][2] + zij;
    }

  } else { // orthogonal box
    for (int id = 1; id <= natom; ++id){
      double xij = atpos[id][0] - ref->atpos[id][0];
      double yij = atpos[id][1] - ref->atpos[id][1];
      double zij = atpos[id][2] - ref->atpos[id][2];

      while (xij >= hx) xij -= lx;
      while (xij < -hx) xij += lx;
    
      while (yij >= hy) yij -= ly;
      while (yij < -hy) yij += ly;
    
      while (zij >= hz) zij -= lz;
      while (zij < -hz) zij += lz;

      atpos[id][0] = ref->atpos[id][0] + xij;
      atpos[id][1] = ref->atpos[id][1] + yij;
      atpos[id][2] = ref->atpos[id][2] + zij;
    }
  }

return;
}

/* -----------------------------------------------------------------------------
 * Private method to apply PBC on a cartesian vector
 * ---------------------------------------------------------------------------*/
void DumpAtom::apply_pbc(double *rij)
{
  if (triclinic) {  // non-orthogonal box
    while (rij[2] >= hz){
      rij[0] -= xz;
      rij[1] -= yz;
      rij[2] -= lz;
    }
    while (rij[2] < -hz){
      rij[0] += xz;
      rij[1] += yz;
      rij[2] += lz;
    }

    while (rij[1] >= hy){
      rij[0] -= xy;
      rij[1] -= ly;
    }
    while (rij[1] < -hy){
      rij[0] += xy;
      rij[1] += ly;
    }

    while (rij[0] >= hx) rij[0] -= lx;
    while (rij[0] < -hx) rij[0] += lx;

  } else { // orthogonal box

    while (rij[0] >= hx) rij[0] -= lx;
    while (rij[0] < -hx) rij[0] += lx;

    while (rij[1] >= hy) rij[1] -= ly;
    while (rij[1] < -hy) rij[1] += ly;

    while (rij[2] >= hz) rij[2] -= lz;
    while (rij[2] < -hz) rij[2] += lz;
  }
return;
}
/* -------------------------------------------------------------------------- */
