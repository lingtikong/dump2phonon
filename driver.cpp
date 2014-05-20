#include "atom.h"
#include "driver.h"
#include "fix_phonon.h"
#include <algorithm>
#include "version.h"

/* --------------------------------------------------------------------------- *
 * Constructor, main driver
 * --------------------------------------------------------------------------- */
Driver::Driver(int narg, char **arg)
{
  first = 1;
  one = ref = NULL;
  flist.clear();
  fdump = NULL;
  char *ctrl; ctrl = NULL;

  // analyse command line options
  int iarg = 1;
  while (narg > iarg){
    if (strcmp(arg[iarg],"-h") == 0){ // help
      help();

    } else { // control file and dump file(s)
      if (ctrl == NULL){
        ctrl = new char [strlen(arg[iarg])+1];
        strcpy(ctrl, arg[iarg]);
      } else {
        fname.assign(arg[iarg]);
        if (find(flist.begin(), flist.end(), fname) == flist.end()) flist.push_back(fname);
      }
    }

    ++iarg;
  }

  // check command line options
  if (ctrl == NULL){
    printf("\nError: insufficient command line options!\n");
    help();
  }

  int ndump = flist.size();
  if (ndump < 1){
    printf("\nError: no dump file passed!\n");
    help();
  }

  // initialize fix-phonon
  FixPhonon *phonon = new FixPhonon(ctrl);
  delete []ctrl;
  if (phonon->status > 0) return;

  // read first/reference frame from dump file(s)
  fname = flist.front();
  FILE *fp = fopen(fname.c_str(), "r");
  if (fp == NULL){
    printf("\nError: file %s not found!\n", fname.c_str());
    help();
  }
  ref = new DumpAtom(fp, fname.c_str());
  fclose(fp);

  if (ref->initialized == 0){
    printf("\nError while reading the first from from dump file %s!\n", fname.c_str());
    help();
  }

  // initialize phonon
  phonon->one = ref;
  if (phonon->init()){
    printf("\nError while initializing the phonon calculations!\n");
    help();
  }
  phonon->setup();

  // now to proceed the real computation
  for (int i = 0; i <= phonon->nskip; ++i) readdump();
  while (one->initialized){
    one->wrap2ref(ref);
    phonon->one = one;

    phonon->end_of_step();
    delete one; one = NULL;

    readdump();
  }

  // finalize
  phonon->post_run();

  fname.clear();
return;
}

/* --------------------------------------------------------------------------- *
 * free memories
 * --------------------------------------------------------------------------- */
Driver::~Driver()
{
  if (ref) delete ref;
  if (one) delete one;
  if (fdump) fclose(fdump);

  fname.clear();
return;
}

/* --------------------------------------------------------------------------- *
 * help
 * --------------------------------------------------------------------------- */
void Driver::help()
{
  printf("\nd2p  version 0.%d, compiled on %s %s\n", VERSION, __DATE__, __TIME__);
  printf("\nPhonons from MD trajectories.\n");
  printf("\nUsage:\n   d2p control-file  dump-file [dump-file2]\n");
  printf("\n");
  exit(0);
}

/* --------------------------------------------------------------------------- *
 * Read one image and assign to "one"
 * --------------------------------------------------------------------------- */
void Driver::readdump()
{
  if (first == 1){
    while( fdump == NULL && flist.size() > 0){
      fname = flist.front(); flist.pop_front();
      fdump = fopen(fname.c_str(), "r");
      if (fdump == NULL){
        printf("\nWarning: faild to open file: %s\n", fname.c_str());
      } else {
        printf("Now to process file %s, takes time...\n", fname.c_str());
      }
    }

    first = 0;
  }

  one = new DumpAtom(fdump, fname.c_str());

  if (one->initialized == 0 && flist.size() > 0){
    while (flist.size() > 0){
      fclose(fdump);
      fname = flist.front(); flist.pop_front();
      fdump = fopen(fname.c_str(), "r");
      if (fdump == NULL){
        printf("\nWarning: faild to open file: %s\n", fname.c_str());
      } else {
        printf("Now to process file %s, takes time...\n", fname.c_str());
        break;
      }
    }

    if (fdump){
      delete one; one = NULL;
      readdump();
    }
  }

return;
}
    
/* --------------------------------------------------------------------------- */
