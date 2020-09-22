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
This is a electrical double layer force model to be used in the particle-bubble interaction.
Contributing author: Linhan Ge (University of Newcastle)
Reference:
Eq.13 from
Yoon, R.H., Flinn, D.H., Rabinovich, Y.I., 1997. Hydrophobic interactions between dissimilar surfaces. J Colloid Interf Sci 185, 363-370. 
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_edl.h"
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

PairEdl::PairEdl(LAMMPS *lmp) : Pair(lmp) {}

/* -----------------------------------------------------------------------*/

PairEdl::~PairEdl()
{
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);

		memory->destroy(cut);
		memory->destroy(epsilon);
		memory->destroy(psi1);
		memory->destroy(psi2);
		memory->destroy(lowcut);
	}
}

/* -----------------------------------------------------------------------*/

void PairEdl::compute(int eflag, int vflag)
{
	int i,j,ii,jj,inum,jnum,itype,jtype;
	double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
	double rsq,r,rinv,V_edl,term1,H,Hinv;
	double psi1ij,psi2ij,epsilonij,cutij,lowcutij;
	double radi,radj,radsum,radtimes;
	double fx,fy,fz;
	int *ilist, *jlist, *numneigh, **firstneigh;
		
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
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx*delx + dely*dely + delz*delz;
			radj = radius[j];
			radsum = radi + radj;
			radtimes = radi*radj;
			jtype = type[j];
			r = sqrt(rsq);
			cutij = cut[itype][jtype];
			lowcutij = lowcut[itype][jtype];
			// H = (r -radsum) > lowcutij ? (r-radsum) : lowcutij;
			H = MAX(r-radsum,lowcutij);
			Hinv = 1.0/H;
			rinv = 1.0/r;
			term1 = radtimes/radsum;
			epsilonij = epsilon[itype][jtype];
			psi1ij = psi1[itype][jtype];
			psi2ij = psi2[itype][jtype];
			if (rsq > radsum*radsum && rsq <= (radsum+cutij)*(radsum+cutij)) {     // do not use cusq
				V_edl = 0.25*epsilonij*term1*(psi1ij*psi1ij+psi2ij*psi2ij)* \
					(
						2*psi1ij*psi2ij/(psi1ij*psi1ij+psi2ij*psi2ij)* \
						log((1+exp(-H/kappainv))/(1-exp(-H/kappainv)))+ \
						log(1-exp(-2*H/kappainv))
					);
				fpair = V_edl * Hinv;
				fx = fpair*delx*rinv;
				fy = fpair*dely*rinv;
				fz = fpair*delz*rinv;
				f[i][0] += fx;
				f[i][1] += fy;
				f[i][2] += fz;
				if (newton_pair || j < nlocal) {
					f[j][0] -= fx;
					f[j][1] -= fy;
					f[j][2] -= fz;
				}
				// set j = nlocal so that only I gets tallied
				if (evflag) ev_tally_xyz(i,nlocal,nlocal,0,0.0,0.0,-fx,-fy,-fz,delx,dely,delz,1,1);
			}
			
		}
	}

	if (vflag_fdotr) virial_fdotr_compute();

}
/* --------------------------------------------------------------------------------

allocate all arrays

---------------------------------------------------------------------------------*/

void PairEdl::allocate()
{
	allocated = 1;
	int n = atom->ntypes;

	memory->create(setflag,n+1,n+1,"pair:setflag");
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			setflag[i][j] = 0;

	memory->create(cutsq,n+1,n+1,"pair:cutsq");
	memory->create(cut,n+1,n+1,"pair:cut");
	memory->create(epsilon,n+1,n+1,"pair:epsilon");
	memory->create(psi1,n+1,n+1,"pair:psi1");
	memory->create(psi2,n+1,n+1,"pair:psi2");
	memory->create(lowcut,n+1,n+1,"pair:lowcut");
}
/* -------------------------------------------------------------------------

global settings

--------------------------------------------------------------------------- */

void PairEdl::settings(int narg, char **arg)
{

	if (narg != 2) error->all(FLERR,"Illegal pair_style command");

	kappainv = force->numeric(FLERR,arg[0]);
	cut_global = force->numeric(FLERR,arg[1]);

	// reset cutoffs that have been explicitly set

	if (allocated) {
		int i,j;
		for (i = 1; i <= atom->ntypes; i++)
			for (j = i+1; j <= atom->ntypes; j++)
				if (setflag[i][j]) cut[i][j] = cut_global;
	}
}
/* ------------------------------------------------------------------------

set coeffs for one or more type pairs

---------------------------------------------------------------------------*/

void PairEdl::coeff(int narg, char **arg)
{
	if (narg < 6 || narg > 7) error->all(FLERR,"Incorrect args for pair coefficients");
	if (!allocated) allocate();

	int ilo,ihi,jlo,jhi;
	force->bounds(arg[0],atom->ntypes,ilo,ihi);
	force->bounds(arg[1],atom->ntypes,jlo,jhi);

	double epsilona_one = force->numeric(FLERR,arg[2]);
	double psi1_one = force->numeric(FLERR,arg[3]);
	double psi2_one = force->numeric(FLERR,arg[4]);
	double lowcut_one = force->numeric(FLERR,arg[5]);

	double cut_one = cut_global;
	if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo,i); j <= jhi; j++) {
			epsilon[i][j] = epsilona_one;
			psi1[i][j] = psi1_one;
			psi2[i][j] = psi2_one;
			lowcut[i][j] = lowcut_one;
			cut[i][j] = cut_one;
			setflag[i][j] = 1;
			count++;
		}
	}

	if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/*------------------------------------------------------------------------

init for one type pair i,j and corresponding j,i

-------------------------------------------------------------------------- */

double PairEdl::init_one(int i, int j)
{
	if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

	epsilon[i][j] = epsilon[i][j];
	psi1[i][j] = psi1[j][i];
	psi2[i][j] = psi2[j][i];
	lowcut[i][j] = lowcut[j][i];

	return cut[i][j];
}
/*------------------------------------------------------------------------

proc 0 writes to restart file

-------------------------------------------------------------------------- */

void PairEdl::write_restart(FILE *fp)
{
	write_restart_settings(fp);

	int i,j;
	for (i = 1; i <= atom->ntypes; i++)
		for (j = i; j <= atom->ntypes; j++) {
			fwrite(&setflag[i][j],sizeof(int),1,fp);
			if (setflag[i][j]) {
				fwrite(&epsilon[i][j],sizeof(double),1,fp);
				fwrite(&psi1[i][j],sizeof(double),1,fp);
				fwrite(&psi2[i][j],sizeof(double),1,fp);
				fwrite(&lowcut[i][j],sizeof(double),1,fp);
				fwrite(&cut[i][j],sizeof(double),1,fp);
			}
		}
}

/*------------------------------------------------------------------------

proc 0 reads from restart file, bcasts

-------------------------------------------------------------------------- */

void PairEdl::read_restart(FILE *fp)
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
					fread(&epsilon[i][j],sizeof(double),1,fp);
					fread(&psi1[i][j],sizeof(double),1,fp);
					fread(&psi2[i][j],sizeof(double),1,fp);
					fread(&lowcut[i][j],sizeof(double),1,fp);
					fread(&cut[i][j],sizeof(double),1,fp);
				}
				MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&psi1[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&psi2[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&lowcut[i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
			}
		}
}

/*------------------------------------------------------------------------
proc 0 writes to restart file
-------------------------------------------------------------------------- */

void PairEdl::write_restart_settings(FILE *fp)
{
	fwrite(&kappainv,sizeof(double),1,fp);
	fwrite(&cut_global,sizeof(double),1,fp);
	//fwrite(&offset_flag,sizeof(int),1,fp);
	fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ------------------------------------------------------------------------
proc 0 reads from restart file, bcasts
------------------------------------------------------------------------ */

void PairEdl::read_restart_settings(FILE *fp)
{
	if (comm->me == 0) {
		fread(&kappainv,sizeof(double),1,fp);
		fread(&cut_global,sizeof(double),1,fp);
		fread(&mix_flag,sizeof(int),1,fp);
	}
	MPI_Bcast(&kappainv,1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
	MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}


/*--------------------------------------------------------- */

double PairEdl::single(int i, int j, int itype,int jtype,
	double rsq,
	double factor_coul, double factor_lj,
	double &fforce)
{
	double r,radi,radj,term1,radsum,radtimes,H,Hinv,rinv;
	double psi1ij,psi2ij,epsilonij,lowcutij;
	double V_edl,fpair;

	double *radius = atom->radius;
	radi = radius[i];
	radj = radius[j];
	radsum = radi + radj;
	radtimes = radi*radj;

	term1 = M_PI*radtimes/radsum;
	epsilonij = epsilon[itype][jtype];
	psi1ij = psi1[itype][jtype];
	psi2ij = psi2[itype][jtype];
	r = sqrt(rsq);
	lowcutij = lowcut[itype][jtype];
	// H = (r -radsum) > lowcutij ? (r-radsum) : lowcutij;
	H = MAX(r-radsum,lowcutij);
	Hinv = 1.0/H;
	rinv = 1/r;
	V_edl = 0.25*epsilonij*term1*(psi1ij*psi1ij+psi2ij*psi2ij)* \
			(
				2*psi1ij*psi2ij/(psi1ij*psi1ij+psi2ij*psi2ij)* \
				log((1+exp(-H/kappainv))/(1-exp(-H/kappainv)))+ \
				log(1-exp(-2*H/kappainv))
			);
	fpair = V_edl * Hinv;
	fforce = fpair*V_edl*rinv;
	return factor_lj*V_edl;
} 

