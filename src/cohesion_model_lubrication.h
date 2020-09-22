/* ----------------------------------------------------------------------
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

	Linhan Ge (University of Newcastle)

	Copyright 2012-     DCS Computing GmbH, Linz
	Copyright 2009-2012 JKU Linz

	This is a lubrication model using Taylor equation.
	Contributor: Linhan Ge
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_LUBRICATION,lubrication,4)
#else

#ifndef COHESION_MODEL_LUBRICATION_H_
#define COHESION_MODEL_LUBRICATION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <math.h>
#include "global_properties.h"
#include "neighbor.h"

namespace MODEL_PARAMS
{
	inline static ScalarProperty* createFluidViscosity(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
	  return fluidViscosityScalar;
	}
	
	inline static MatrixProperty* createMinSeparationDist(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  MatrixProperty* minSeparationDistMatrix = MODEL_PARAMS::createPerTypePairProperty(registry, "minSeparationDist", caller);
	  return minSeparationDistMatrix;
	}
    
	inline static ScalarProperty* createMaxSeparationDistRatio(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* maxSeparationDistRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistRatio", caller);
	  return maxSeparationDistRatioScalar;
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
	  CohesionModelBase(lmp, hsetup, cmb), fluidViscosity(0.0), maxSeparationDistRatio(0.), minSeparationDist(NULL)
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
		registry.registerProperty("minSeparationDist", &MODEL_PARAMS::createMinSeparationDist);
		registry.registerProperty("maxSeparationDistRatio", &MODEL_PARAMS::createMaxSeparationDistRatio);
		registry.connect("fluidViscosity", fluidViscosity,"cohesion_model lubrication");
		registry.connect("minSeparationDist", minSeparationDist,"cohesion_model lubrication");
		registry.connect("maxSeparationDistRatio", maxSeparationDistRatio,"cohesion_model lubrication");

		// error checks on coarsegraining
		if(force->cg_active())
			error->cg(FLERR,"cohesion model lubrication");

		neighbor->register_contact_dist_factor(maxSeparationDistRatio); 
		if(maxSeparationDistRatio < 1.0)
			error->one(FLERR,"\n\ncohesion model lubrication requires maxSeparationDistanceRatio >= 1.0. Please increase this value.\n");
	}

	inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
	void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
	void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

	void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
	   if(sidata.contact_flags) *sidata.contact_flags &= ~CONTACT_COHESION_MODEL;
	}

	void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
	{
	   
	   if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;

	   scdata.has_force_update = true;

	   const int itype = scdata.itype;
       const int jtype = scdata.jtype;

	   if (!scdata.is_wall) {

		  const double rsq = scdata.rsq;
		  const double r = sqrt(rsq);
		  const double rinv =  1.0/r;
		  const double radsum = scdata.radsum;
		  const double radi = scdata.radi;
		  const double radj = scdata.radj;
		  const double rEff = radi*radj / radsum;
		  double d = r - radsum;
		  
		  d = d > minSeparationDist[itype][jtype] ? d : minSeparationDist[itype][jtype];
			  
		  const double dx = scdata.delta[0];
		  const double dy = scdata.delta[1];
		  const double dz = scdata.delta[2];
		  const double enx = dx * rinv;
		  const double eny = dy * rinv;
		  const double enz = dz * rinv;
		  // relative translational velocity
		  const double vr1 = scdata.v_i[0] - scdata.v_j[0];
		  const double vr2 = scdata.v_i[1] - scdata.v_j[1];
		  const double vr3 = scdata.v_i[2] - scdata.v_j[2];
		  const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
		  
		  const double F_lubrication = -6*M_PI*fluidViscosity*vn*rEff*rEff/d;
			  
		  const double fx = F_lubrication * enx;       //en represent the normal direction vector, en[0] is the x coordinate
		  const double fy = F_lubrication * eny;				 
		  const double fz = F_lubrication * enz;				 

		  i_forces.delta_F[0] += fx;
		  i_forces.delta_F[1] += fy;
		  i_forces.delta_F[2] += fz;

		  j_forces.delta_F[0] -= fx;
		  j_forces.delta_F[1] -= fy;
		  j_forces.delta_F[2] -= fz;
	  }

	  if(scdata.is_wall) {
		
		double d = scdata.nonConDeltan;                             // deltan is the distance to the wall if scdata.wall = true
		d = d > minSeparationDist[itype][jtype] ? d : minSeparationDist[itype][jtype];
		const double rinv =  1.0/scdata.nonConr;
		const double enx = scdata.delta[0] * rinv;
		const double eny = scdata.delta[1] * rinv;
		const double enz = scdata.delta[2] * rinv;
		const double rEff =  scdata.radi;
						
		const double vr1 = scdata.v_i[0] - scdata.v_j[0];
		const double vr2 = scdata.v_i[1] - scdata.v_j[1];
		const double vr3 = scdata.v_i[2] - scdata.v_j[2];

		const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
		
		const double F_lubrication = -6*M_PI*fluidViscosity*vn*rEff*rEff/d;
			
		const double fx = F_lubrication * enx;    			//en represent the normal direction vector, en[0] is the x coordinate
		const double fy = F_lubrication * eny;				 
		const double fz = F_lubrication * enz;				 
			
		i_forces.delta_F[0] += fx;
		i_forces.delta_F[1] += fy;
		i_forces.delta_F[2] += fz;
	  }
	}
  
  private:
	double fluidViscosity, maxSeparationDistRatio;
	double ** minSeparationDist;
	bool tangentialReduce_;
  };
 }
}

#endif 
#endif