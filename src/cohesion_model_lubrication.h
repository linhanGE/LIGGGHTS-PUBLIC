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
COHESION_MODEL(COHESION_LUBRICATION,lubrication,4)
#else

#ifndef COHESION_MODEL_LUBRICATION_H_
#define COHESION_MODEL_LUBRICATION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <math.h>

namespace MODEL_PARAMS
{
	inline static ScalarProperty* createFluidViscosity(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
	  return fluidViscosityScalar;
	}
	
	inline static ScalarProperty* createCutoffDistance(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* cutoffDistanceScalar = MODEL_PARAMS::createScalarProperty(registry, "cutoffDistance", caller);
	  return cutoffDistanceScalar;
	}
}

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_LUBRICATION> : public CohesionModelBase {
  public:
	static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;
	
	CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
	  CohesionModelBase(lmp, hsetup, cmb), fluidViscosity(0.0), cutoffDistance(0.0)
	{
		
	}

	void registerSettings(Settings & settings)
	{
	  settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
	}
	
	inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

	void connectToProperties(PropertyRegistry & registry)
	{
		registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosity);
		registry.registerProperty("cutoffDistance", &MODEL_PARAMS::createCutoffDistance);
		registry.connect("fluidViscosity", fluidViscosity,"cohesion_model lubrication");
		registry.connect("cutoffDistance", cutoffDistance,"cohesion_model lubrication");
		// error checks on coarsegraining
		if(force->cg_active())
			error->cg(FLERR,"cohesion model lubrication");
	}

	void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
	  const double r = sidata.r;
	  const double rsq = sidata.rsq;
	  //const int itype = sidata.itype;
	  //const int jtype = sidata.jtype;
	  //const double rhoi = sidata.densityi;
	  //const double rhoj = sidata.densityj;
	  const double radsum = sidata.radsum;
	  const double radi = sidata.radi;
	  const double radj = sidata.is_wall ? radi : sidata.radj;
	  const double rEff = sidata.is_wall ? radi : radi*radj / (radi+radj);
	  double d = 0;

	  if ( rsq > radsum * radsum) {
	  if ( r - radsum < cutoffDistance) d = cutoffDistance;
	  else {
		  d = r -radsum;	
	  }
	   
	  const double Fn = -6*M_PI*fluidViscosity*sidata.vn*rEff*rEff/d;   
	  
	  if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

	  // apply normal force
	  if(sidata.is_wall) {
		const double Fn_ = Fn;
		i_forces.delta_F[0] += Fn_ * sidata.en[0];
		i_forces.delta_F[1] += Fn_ * sidata.en[1];
		i_forces.delta_F[2] += Fn_ * sidata.en[2];
	  } else {
		const double fx = Fn * sidata.en[0];    			//en represent the normal direction vector, en[0] is the x coordinate
		const double fy = Fn * sidata.en[1];				 
		const double fz = Fn * sidata.en[2];				 

		i_forces.delta_F[0] += fx;
		i_forces.delta_F[1] += fy;
		i_forces.delta_F[2] += fz;

		j_forces.delta_F[0] -= fx;
		j_forces.delta_F[1] -= fy;
		j_forces.delta_F[2] -= fz;
	  }
	}
	}

	inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
	void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
	void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

	void surfacesClose(SurfacesCloseData& scdata, ForceData&, ForceData&)
	{
		if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;

	}
  
  private:
	double fluidViscosity,cutoffDistance;
	bool tangentialReduce_;
  };
}
}

#endif 
#endif