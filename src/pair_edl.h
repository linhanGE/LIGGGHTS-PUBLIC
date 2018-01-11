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

/* ----------------------------------------------------------------------
This is an electrical double layer force modle to be used in the particle-bubble interaction.
contributor: Linhan Ge (University of Newcastle)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(edl,PairEdl)

#else

#ifndef LMP_PAIR_EDL_H
#define LMP_PAIR_EDL_H

#include "pair.h"

namespace LAMMPS_NS {

	class PairEdl : public Pair {
	public:
		PairEdl(class LAMMPS *);
		virtual ~PairEdl();
		virtual void compute(int, int);
		void settings(int, char **);
		void coeff(int, char **);
		virtual double init_one(int, int);
		void write_restart(FILE *);
		void read_restart(FILE *);
		void write_restart_settings(FILE *);
		void read_restart_settings(FILE *);
		virtual double single(int, int, int, int, double, double, double, double &);

	protected:
		double cut_global;
		double kappa;
		double **cut;
		double **lowcut,**epsilona,**psi1,**psi2;

		void allocate();

	};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
