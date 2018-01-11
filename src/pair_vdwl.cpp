/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   contributor: Linhan Ge (University of Newcastle)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "pair_vdwl.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

PairVdwl::PairVdwl(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairVdwl::~PairVdwl()
{
  if (allocated) {
	memory->destroy(setflag);
	memory->destroy(cutsq);
	
	memory->destroy(A132);
	memory->destroy(lowcut);
	memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairVdwl::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,force_vdwl;
  double radi,radj,radsum,radtimes,r,factor_lj;
  double A132ij,H,H2,lowcutij;
  double term1;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
	i = ilist[ii];
	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];
	radi = radius[i];
	itype = type[i];
	jlist = firstneigh[i];
	jnum = numneigh[i];

	for (jj = 0; jj < jnum; jj++) {
	  j = jlist[jj];
	  factor_lj = special_lj[sbmask(j)];
	  j &= NEIGHMASK;

	  delx = xtmp - x[j][0];
	  dely = ytmp - x[j][1];
	  delz = ztmp - x[j][2];
	  rsq = delx*delx + dely*dely + delz*delz;
	  r = sqrt(rsq);
	  radj = radius[j];
	  radsum = radi + radj;
	  jtype = type[j];
	  lowcutij = lowcut[itype][jtype];
	  H = (r -radsum) > lowcutij ? (r-radsum) : lowcutij;
	  H2 = H*H;
	  double Hinv = 1.0/H;
	  double rinv = 1.0/r;
	  if (rsq >= radsum*radsum) {
		A132ij = A132[itype][jtype];
		term1 = radtimes/radsum;
		force_vdwl = -A132ij/(6*H)*term1;
		fpair = factor_lj*force_vdwl*Hinv*rinv;
		
		f[i][0] += delx*fpair;
		f[i][1] += dely*fpair;
		f[i][2] += delz*fpair;
		if (newton_pair || j < nlocal) {
		  f[j][0] -= delx*fpair;
		  f[j][1] -= dely*fpair;
		  f[j][2] -= delz*fpair;
		}

		if (eflag) {
		  evdwl = -A132ij/(6*H)*term1;
		}

		if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
							 evdwl,0.0,f[i][0],f[i][1],f[i][2],delx,dely,delz);
	  }
	}
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairVdwl::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
	for (int j = i; j <= n; j++)
	  setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(A132,n+1,n+1,"pair:A132");
  memory->create(lowcut,n+1,n+1,"pair:lowcut");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairVdwl::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
	int i,j;
	for (i = 1; i <= atom->ntypes; i++)
	  for (j = i+1; j <= atom->ntypes; j++)
		if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairVdwl::coeff(int narg, char **arg)
{
  if (narg != 4 || narg != 5)
	error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double A132_one = force->numeric(FLERR,arg[2]);
  double lowcut_one = force->numeric(FLERR,arg[3]);
  
  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
	for (int j = MAX(jlo,i); j <= jhi; j++) {
	  A132[i][j] = A132_one;
	  cut[i][j] = cut_one;
	  lowcut[i][j] = lowcut_one;
	  setflag[i][j] = 1;
	  count++;
	}
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairVdwl::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  A132[j][i] = A132[i][j];
  lowcut[j][i] = lowcut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairVdwl::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
	for (j = i; j <= atom->ntypes; j++) {
	  fwrite(&setflag[i][j],sizeof(int),1,fp);
	  if (setflag[i][j]) {
		fwrite(&A132[i][j],sizeof(double),1,fp);
		fwrite(&lowcut[i][j],sizeof(double),1,fp);
		fwrite(&cut[i][j],sizeof(double),1,fp);
	  }
	}
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairVdwl::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
	for (j = i; j <= atom->ntypes; j++) {
	  if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
	  MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
	  if (setflag[i][j]) {
		if (me == 0) {
		  fread(&A132[i][j],sizeof(double),1,fp);
		  fread(&lowcut[i][j],sizeof(double),1,fp);
		  fread(&cut[i][j],sizeof(double),1,fp);
		}
		MPI_Bcast(&A132[i][j],1,MPI_DOUBLE,0,world);
		MPI_Bcast(&lowcut[i][j],1,MPI_DOUBLE,0,world);
		MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
	  }
	}
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairVdwl::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairVdwl::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
	fread(&cut_global,sizeof(double),1,fp);
	fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairVdwl::single(int i, int j, int itype, int jtype, double rsq,
						  double factor_coul, double factor_lj,
						  double &fforce)
{
  double *radius =atom->radius;

  double radi = radius[i];
  double radj = radius[j];
  double radsum = radi + radj;
  double radtimes = radi*radj;
  
  double A132ij = A132[itype][jtype];
  double lowcutij = lowcut[itype][jtype];
  double r = sqrt(rsq);
  double term1 = radtimes/radsum;
  double H = (r -radsum) > lowcutij ? (r-radsum) : lowcutij;
  double Hinv = 1.0/H;
  double rinv = 1.0/r;
  double force_vdwl = -A132ij/(6*H)*term1;
  fforce = factor_lj*force_vdwl*Hinv*rinv;
  double phi_vdwl = -A132ij/(6*H)*term1;
  
  return factor_lj*phi_vdwl;
}