/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_evaporate.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "force.h"
#include "group.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEvaporate::FixEvaporate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix evaporate command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  // nevery = atoi(arg[3]);
  nflux = atoi(arg[3]);
  char * seed = arg[4];
  keepGroup = group->find(arg[5]);
  
  int iarg = 6;
  if (strcmp(arg[iarg],"region") == 0)
  {
    if (iarg+2 > narg) error->all(FLERR,"Illegal fix evaporate command");
    iregion = domain->find_region(arg[iarg+1]);
    int n = strlen(arg[iarg+1]) + 1;
    idregion = new char[n];
    strcpy(idregion,arg[iarg+1]);
  } else error->all(FLERR,"Illegal fix evaporate command");
  
  if (nflux <= 0)
    error->all(FLERR,"Illegal fix evaporate command");

  if (seed <= 0) error->all(FLERR,"Illegal fix evaporate command");

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);
 
  // set up reneighboring
  /*force_reneighbor = 1;
  next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;*/
  ndeleted = 0;

  nmax = 0;
  list = NULL;
  mark = NULL;
}

/* ---------------------------------------------------------------------- */

FixEvaporate::~FixEvaporate()
{
  delete [] idregion;
  delete random;
  memory->destroy(list);
  memory->destroy(mark);
}

/* ---------------------------------------------------------------------- */

int FixEvaporate::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEvaporate::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);

  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot evaporate atoms in atom_modify first group");
  }
}

/* ----------------------------------------------------------------------
   perform particle deletion
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixEvaporate::pre_exchange()
{
  int i,iwhichglobal,iwhichlocal;
  int ndel,ndeltopo[4];

  // if (update->ntimestep != next_reneighbor) return;

  // grow list and mark arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(list);
    memory->destroy(mark);
    nmax = atom->nmax;
    memory->create(list,nmax,"evaporate:list");
    memory->create(mark,nmax,"evaporate:mark");
  }

  // ncount = # of deletable atoms in region that I own
  // nall = # on all procs
  // nbefore = # on procs before me
  // list[ncount] = list of local indices of atoms I can delete

  double **x = atom->x;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ncount = 0;
  int keepGroupbit = group->bitmask[keepGroup];

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & keepGroupbit)
      continue;
    if (mask[i] & groupbit)
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
        list[ncount++] = i;
  }

  int nall,nbefore;
  MPI_Allreduce(&ncount,&nall,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ncount,&nbefore,1,MPI_INT,MPI_SUM,world);
  nbefore -= ncount;

  // ndel = total # of atom deletions, in or out of region
  // ndeltopo[1,2,3,4] = ditto for bonds, angles, dihedrals, impropers
  // mark[] = 1 if deleted

  ndel = 0;
  for (i = 0; i < nlocal; i++) mark[i] = 0;

  // atomic deletions
  // choose atoms randomly across all procs and mark them for deletion
  // shrink eligible list as my atoms get marked
  // keep ndel,ncount,nall,nbefore current after each atom deletion

  while (nall && ndel < nflux) {
    iwhichglobal = static_cast<int> (nall*random->uniform());
    if (iwhichglobal < nbefore) nbefore--;
    else if (iwhichglobal < nbefore + ncount) {
      iwhichlocal = iwhichglobal - nbefore;
      mark[list[iwhichlocal]] = 1;
      list[iwhichlocal] = list[ncount-1];
      ncount--;
    }
    ndel++;
    nall--;
  }
  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms

  AtomVec *avec = atom->avec;

  for (i = nlocal-1; i >= 0; i--) {
    if (mark[i]) {
      avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
    }
  }

  // compress id
  /*for (i = 0; i < nlocal; i++) tag[i] = 0;
  atom->tag_extend();*/  

  // reset global natoms and bonds, angles, etc
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  atom->natoms -= ndel;

  if (ndel && atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // statistics

  ndeleted += ndel;
  // next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
   return number of deleted particles
------------------------------------------------------------------------- */

double FixEvaporate::compute_scalar()
{
  return 1.0*ndeleted;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixEvaporate::memory_usage()
{
  double bytes = 2*nmax * sizeof(int);
  return bytes;
}
