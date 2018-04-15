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
	
	inline static ScalarProperty* createMinSeparationDist(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* minSeparationDistScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDist", caller);
	  return minSeparationDistScalar;
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
	  CohesionModelBase(lmp, hsetup, cmb), 
      coeffMu(NULL),
      fluidViscosity(0.),
      minSeparationDist(0.),
      maxSeparationDistRatio(0.),
      particleOnly(true)
	{
		
	}

	inline void registerSettings(Settings& settings) 
    {
        settings.registerOnOff("particleOnly", particleOnly, true);
    }
    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}
    
	void connectToProperties(PropertyRegistry & registry)
	{
		registry.registerProperty("coeffMu", &MODEL_PARAMS::createCoeffMu);
		registry.registerProperty("minSeparationDist", &MODEL_PARAMS::createMinSeparationDist);
		registry.registerProperty("maxSeparationDistRatio", &MODEL_PARAMS::createMaxSeparationDistRatio);
		registry.connect("coeffMu", coeffMu,"cohesion_model lubrication");
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
       const int i = scdata.i;
	   const int j = scdata.j;
       const int itype = scdata.itype;
       const int jtype = scdata.jtype;
	   const double radi = scdata.radi;
	   const double radj = scdata.is_wall ? radi : scdata.radj;
	   const double r = sqrt(scdata.rsq);
	   const double radsum = scdata.radsum;
	   const double dist = scdata.is_wall ? r - radi : r - (radi + radj);
	   const double rEff = scdata.is_wall ? radi : radi*radj / radsum;
	   const double d = dist > minSeparationDist ? dist : minSeparationDist;

       const double fluidViscosity = coeffMu[itype][jtype];

	   double **v = atom->v;
	   // calculate vn and vt since not in struct
	   const double rinv = 1.0 / r;
	   const double dx = scdata.delta[0];
	   const double dy = scdata.delta[1];
	   const double dz = scdata.delta[2];
	   const double enx = dx * rinv;
	   const double eny = dy * rinv;
	   const double enz = dz * rinv;

	   // relative translational velocity
	   const double vr1 = v[i][0] - v[j][0];
	   const double vr2 = v[i][1] - v[j][1];
	   const double vr3 = v[i][2] - v[j][2];

	   // normal component
	   const double vn = vr1 * enx + vr2 * eny + vr3 * enz;

	   double F_lubrication = -6*M_PI*fluidViscosity*vn*rEff*rEff/d;

       if (scdata.is_wall && particleOnly) F_lubrication = 0.;
			  
	   const double fx = F_lubrication * enx;    			//en represent the normal direction vector, en[0] is the x coordinate
	   const double fy = F_lubrication * eny;				 
	   const double fz = F_lubrication * enz;				 

	   i_forces.delta_F[0] += fx;
	   i_forces.delta_F[1] += fy;
	   i_forces.delta_F[2] += fz;
       
       j_forces.delta_F[0] -= fx;
	   j_forces.delta_F[1] -= fy;
	   j_forces.delta_F[2] -= fz;
	}
  
  private:
    double ** coeffMu;
	double fluidViscosity,minSeparationDist, maxSeparationDistRatio;	
    bool   particleOnly;
  };
 }
}

#endif 
#endif