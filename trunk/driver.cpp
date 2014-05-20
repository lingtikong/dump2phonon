#include "atom.h"
#include "driver.h"
#include "fix_phonon.h"
#include <algorithm>

/* --------------------------------------------------------------------------- *
 * Constructor, main driver
 * --------------------------------------------------------------------------- */
Driver::Driver(int narg, char **arg)
{
  first = 1;
  one = ref = NULL;
  flist.clear();
  fdump = NULL;

  // analyze command line options
  if (narg < 3){
    printf("\nError: insufficient command line options!\n");
    help();
  }

  // strip dump file name(s)
  for (int iarg = 2; iarg < narg; ++iarg){
    fname.assign(arg[iarg]);
    if (find(flist.begin(), flist.end(), fname) == flist.end()) flist.push_back(fname);
  }
  int ndump = flist.size();
  if (ndump < 1){
    printf("\nError: no dump file passed!\n");
    help();
  }

  // initialize fix-phonon
  FixPhonon *phonon = new FixPhonon(arg[1]);
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
    delete one;

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
  if (fdump) fclose(fdump);

  fname.clear();
return;
}

/* --------------------------------------------------------------------------- *
 * help
 * --------------------------------------------------------------------------- */
void Driver::help()
{
  printf("\nd2p: phonons from MD trajectories.\n");
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
        break;
      }
    }

    if (fdump){
      delete one;
      readdump();
    }
  }

return;
}
    
/* --------------------------------------------------------------------------- */
