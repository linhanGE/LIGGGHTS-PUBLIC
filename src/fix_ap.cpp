/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "fix_ap.h"

#include "atom.h"
#include "fix_property_atom.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "pair_gran.h"
#include <stdlib.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAp::FixAp(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix fold needs per particle radius and mass");

  if (narg < 5)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;

  fix_ap = NULL;
}

/* ---------------------------------------------------------------------- */

void FixAp::post_create()
{
  // register directional flux
  fix_ap = static_cast<FixPropertyAtom*>(modify->find_fix_property("ap","property/atom","vector",3,0,this->style,false));
  if(!fix_ap)
  {
    const char* fixarg[11];
    fixarg[0]="ap";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="ap";
    fixarg[4]="vector";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_ap = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
  }

  fix_ap = static_cast<FixPropertyAtom*>(modify->find_fix_property("ap","property/atom","vector",0,0,style));

  if(!fix_ap)
    error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

void FixAp::updatePtrs()
{
  ap = fix_ap->array_atom;
}

/* ---------------------------------------------------------------------- */

void FixAp::init()
{
  
  if (!atom->radius_flag || !atom->rmass_flag)
    error->fix_error(FLERR,this,"must use a granular atom style ");

    // check if a fix of this style already exists
  if(modify->n_fixes_style(style) > 1)
    error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  fix_ap = static_cast<FixPropertyAtom*>(modify->find_fix_property("ap","property/atom","vector",0,0,style));

  if(!fix_ap)
    error->one(FLERR,"internal error");

  updatePtrs();

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  {

     //if (mask[i] & groupbit)
     {
        ap[i][0] = 0.;
        ap[i][1] = 0.;
        ap[i][2] = 0.;
     }
  }
}

/* ---------------------------------------------------------------------- */

int FixAp::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAp::end_of_step(int vflag)
{
  updatePtrs();

  double **f = atom->f;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  {

     //if (mask[i] & groupbit)
     {
        ap[i][0] = f[i][0]/rmass[i];
        ap[i][1] = f[i][1]/rmass[i];
        ap[i][2] = f[i][2]/rmass[i];
     }
  }
}