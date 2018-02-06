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

	Christoph Kloss (DCS Computing GmbH, Linz)
	Christoph Kloss (JKU Linz)
	Richard Berger (JKU Linz)

	Copyright 2012-     DCS Computing GmbH, Linz
	Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_PB,pb,9)
#else

#ifndef COHESION_MODEL_PB_H_
#define COHESION_MODEL_PB_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_PB> : public CohesionModelBase {
  public:
	CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
		CohesionModelBase(lmp, hsetup, c),
		surfaceTension(0.0),
		contactAngle(NULL)
	{
		
	}

	void registerSettings(Settings& settings) 
	{
		settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
	}

	inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

	void connectToProperties(PropertyRegistry & registry)
	{
		registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
		registry.registerProperty("contactAngle", &MODEL_PARAMS::createContactAngle);
		registry.connect("surfaceTension",surfaceTension ,"cohesion_model pb");
		registry.connect("contactAngle", contactAngle,"cohesion_model pb");
		// error checks on coarsegraining
		if(force->cg_active())
			error->cg(FLERR,"cohesion model pb");
	}

	void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
	  //r is the distance between the sphere's centers
	  const double r = sidata.r;
	  const double ri = sidata.radi;
	  const double rj = sidata.radj;
	  const double deltan = sidata.deltan;
	  const double Ri = 1/(2*r)*sqrt((-r+ri-rj)*(-r-ri+rj)*(-r+ri+rj)*(r+ri+rj));
	  const double rhoi = sidata.densityi;
	  const double rhoj = sidata.densityj;

	  double rb = 0, rp = 0;
	  if (rhoi > 0.1 && rhoj> 0.1) {        // make sure SI unit is used
		  rp = rhoi > 10 ? ri : rj;
		  rb = rhoi > 10 ? rj : ri;		
	  } else {
		  rp = rhoi >= 1 ? ri : rj;
		  rb = rhoi >= 1 ? rj : ri;		
	  }
	  //const double sin_alpha = Ri/rp;
	  const double hp = sqrt(rp * rp - Ri * Ri);
	  const double hb = sqrt(rb * rb - Ri * Ri);
	  
	  const double sp = 2*M_PI*rp*hp;
	  const double sb = 2*M_PI*rb*hb;

	  const double contactEffAngle = rp == rhoi ? contactAngle[sidata.itype] : contactAngle[sidata.itype];

	  const double wa = surfaceTension*(sb-sp*cos(contactEffAngle));

	  const double F_ad = -wa/deltan;
	  if(tangentialReduce_) sidata.Fn += F_ad; 

	  if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
	    
	    const double fx = F_ad * sidata.en[0];
        const double fy = F_ad * sidata.en[1];
        const double fz = F_ad * sidata.en[2];

	  	i_forces.delta_F[0] += fx;
		i_forces.delta_F[1] += fy;
		i_forces.delta_F[2] += fz;

		j_forces.delta_F[0] -= fx;
		j_forces.delta_F[1] -= fy;
		j_forces.delta_F[2] -= fz;
	}

	inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
	void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
	void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

	void surfacesClose(SurfacesCloseData& scdata, ForceData&, ForceData&)
	{
		if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;

	}

  private:
	double surfaceTension;
	double * contactAngle;
	bool tangentialReduce_;
  };
}
}

#endif // COHESION_MODEL_PB_H_
#endif
