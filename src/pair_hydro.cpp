/*----------------------------------------------------------------------------
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov
Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.
See the README file in the top-level LAMMPS directory.
--------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
This is a hydrophobic force model to be used in the particle-bubble interaction.
Contributor: Linhan Ge (University of Newcastle)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_hydro.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "math_special.h"

using namespace LAMMPS_NS;

using namespace MathSpecial;
/*------------------------------------------------------------------------*/

PairHydro::PairHydro(LAMMPS *lmp) : Pair(lmp) {}

/* -----------------------------------------------------------------------*/

PairHydro::~PairHydro()
{
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);
		memory->destroy(lamda);
		memory->destroy(k);
		memory->destroy(lowcut);
		memory->destroy(cut);
	}
}

/* -----------------------------------------------------------------------*/

void PairHydro::compute(int eflag, int vflag)
{
	int i,j,ii,jj,inum,jnum,itype,jtype;
	double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
	double rsq,force_hydro;
	double radi,radj,radsum,radtimes,r,factor_lj,rinv;
	double lamdaij,H,H2,kij,lowcutij;
	double term1;
	int *ilist,*jlist,*numneigh,**firstneigh;

	evdwl =0.0;

	if (eflag || vflag) ev_setup(eflag,vflag);
	else evflag = vflag_fdotr = 0;

	double **x = atom->x;
	double **f = atom->f;
	double *radius = atom->radius;
	int *type = atom->type;
	int nlocal = atom->nlocal;
	int newton_pair = force->newton_pair;

	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	// loop over neighbors of my atoms

	for (ii =0; ii < inum; ii++) {
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
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx*delx + dely*dely + delz*delz;
			radj = radius[j];
			radsum = radi + radj;
			radtimes = radi*radj;
			jtype = type[j];
			r = sqrt(rsq);
			rinv = 1.0/r;
			lowcutij = lowcut[itype] [jtype];
			lamdaij = lamda[itype][jtype];
			kij = k[itype][jtype];
			H = (r -radsum) > lowcutij ? (r-radsum) : lowcutij;  // lower cut-off distance
			H2 = H * H;
			if ( H2 < cutsq[itype][jtype] ) {               // active when overlap
				term1 = radtimes/radsum;                    // harmonic mean of the radius
				force_hydro = -term1*kij*exp(-H/lamdaij);   // pay attention to the sign
				double fx = delx*force_hydro*rinv;
				double fy = dely*force_hydro*rinv;
				double fz = delz*force_hydro*rinv;

				f[i][0] += fx;
				f[i][1] += fy;
				f[i][2] += fz;
				if (newton_pair || j < nlocal) {
					f[j][0] -= fx;
					f[j][1] -= fy;
					f[j][2] -= fz;
				}

				if (eflag) {
					evdwl = -term1*kij*exp(-H/lamdaij)*H;

					if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
						evdwl,0,f[i][0],f[i][1],f[i][2],delx,dely,delz);
				}
			}
		}
		if (vflag_fdotr) virial_fdotr_compute();
	}
}

/* --------------------------------------------------------------------------------

allocate all arrays

---------------------------------------------------------------------------------*/

void PairHydro::allocate()
{
	allocated = 1;
	int n = atom->ntypes;

	memory->create(setflag,n+1,n+1,"pair:setflag");
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			setflag[i][j] = 0;

	memory->create(cutsq,n+1,n+1,"pair:cutsq");
	memory->create(cut,n+1,n+1,"pair:cut");
	memory->create(lamda,n+1,n+1,"pair:lamda");
	memory->create(k,n+1,n+1,"pair:k");
	memory->create(lowcut,n+1,n+1,"pair:lowcut");
}

/* -------------------------------------------------------------------------

global settings

--------------------------------------------------------------------------- */

void PairHydro::settings(int narg, char **arg)
{
	if (narg != 1) error->all(FLERR,"Illegal pair_style command");
	cut_global = atof(arg[0]);

	// reset cutoffs that have been explicitly set

	if (allocated) {
		int i, j;
		for (i = 1; i <= atom->ntypes; i++)
			for (j = i+1; j <= atom->ntypes; j++)
				if (setflag[i][j]) {
					cut[i][j] = cut_global;
				}
	}
}
/* ------------------------------------------------------------------------

set coeffs for one or more type pairs

---------------------------------------------------------------------------*/

void PairHydro::coeff(int narg, char **arg)
{
	if (narg != 5 && narg != 6)
		error->all(FLERR,"Incorrect args for pair coefficients");
	if (!allocated) allocate();

	int ilo,ihi,jlo,jhi;
	force->bounds(arg[0],atom->ntypes,ilo,ihi);	// see details in force.cpp
	force->bounds(arg[1],atom->ntypes,jlo,jhi);

	double lamda_one = atof(arg[2]);	// atof converts string to double
	double k_one = atof(arg[3]);
	double lowcut_one = atof(arg[4]);

	double cut_one = cut_global;	// cut global set in: PairHydro::settings(int narg, char **arg)
	if (narg == 6) cut_one = atof(arg[5]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo,i); j <= jhi; j++) {
			lamda[i][j] = lamda_one;
			k[i][j] = k_one;
			lowcut[i][j] = lowcut_one;
			cut[i][j] = cut_one;
			setflag[i][j] = 1;	// 0/1 = whether each i,j has been set

			count++;

		}
	}

	if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}


/*------------------------------------------------------------------------

init for one type pair i,j and corresponding j,i

-------------------------------------------------------------------------- */

double PairHydro::init_one(int i, int j)
{
	if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

	lamda[j][i] = lamda[i][j];
	k[j][i] = k[i][j];
	lowcut[j][i] = lowcut[i][j];

	return cut[i][j];
}

/*------------------------------------------------------------------------

proc 0 writes to restart file

-------------------------------------------------------------------------- */

void PairHydro::write_restart(FILE *fp)
{
	write_restart_settings(fp);

	int i,j;
	for (i = 1; i <= atom->ntypes; i++)
		for (j = i; j <= atom->ntypes; j++) {
			fwrite(&setflag[i][j],sizeof(int),1,fp);
			if (setflag[i][j]) {
				fwrite(&lamda[i][j],sizeof(double),1,fp);
				fwrite(&k[i][j],sizeof(double),1,fp);
				fwrite(&lowcut[i][j],sizeof(double),1,fp);
				fwrite(&cut[i][j],sizeof(double),1,fp);
			}
		}
}

/*------------------------------------------------------------------------

proc 0 reads from restart file, bcasts

-------------------------------------------------------------------------- */

void PairHydro::read_restart(FILE *fp)
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
					fread(&lamda[i][j],sizeof(double),1,fp);
					fread(&k[i][j],sizeof(double),1,fp);
					fread(&lowcut[i][j],sizeof(double), 1,fp);
					fread(&cut[i][j],sizeof(double),1,fp);
				}
				MPI_Bcast(&lamda[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&k[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&lowcut[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
			}
		}
}

/*------------------------------------------------------------------------
proc 0 writes to restart file
-------------------------------------------------------------------------- */

void PairHydro::write_restart_settings(FILE *fp)
{
	fwrite(&cut_global,sizeof(double),1,fp);
	fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ------------------------------------------------------------------------
proc 0 reads from restart file, bcasts
------------------------------------------------------------------------ */

void PairHydro::read_restart_settings(FILE *fp)
{
	int me = comm->me;
	if (me == 0) {
		fread(&cut_global,sizeof(double),1,fp);
		fread(&mix_flag,sizeof(int),1, fp);
	}
	MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
	MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/*--------------------------------------------------------- */

double PairHydro::single(int i, int j, int itype,int jtype,
	double rsq,
	double factor_coul, double factor_lj,
	double &fforce)
{
	double phi_hydro,r; // separation distance
	double force_hydro;
	double radsum,radtimes;
	double lamdaij,H,kij,lowcutij,radi,radj;
	double term1;

	double *radius = atom->radius;

	radi = radius[i];
	radj = radius[j];
	radsum = radi + radj;
	radtimes = radi*radj;
	lamdaij = lamda[itype][jtype];
	kij = k[itype][jtype];
	lowcutij = lowcut[itype][jtype];
	r = sqrt(rsq);
	term1 = radtimes/radsum;
	H = (r -radsum) > lowcutij ? (r-radsum) : lowcutij; // lower cut-off distance
	double rinv = 1.0/r;
	force_hydro = -term1*kij*exp(-H/lamdaij);
	fforce = factor_lj*force_hydro*rinv;
	phi_hydro = -term1*kij*exp(-H/lamdaij)*H;
	return factor_lj*phi_hydro;
} 

