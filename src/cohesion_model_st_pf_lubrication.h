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

	This cohesion model accounts for the particle-bubble adhesion force and lubrication force.
	Contributor: Linhan Ge
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_ST_PF_LUBRICATION,st/pf/lubrication,9)
#else

#ifndef COHESION_MODEL_ST_PF_LUBRICATION_H_
#define COHESION_MODEL_ST_PF_LUBRICATION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include "global_properties.h"
#include <cmath>
#include <math.h>
#include "neighbor.h"

namespace MODEL_PARAMS
{
	inline static ScalarProperty* createFluidViscosityPb(PropertyRegistry & registry, const char * caller, bool sanity_checks) // Pb stands for particle-bubble
	{
	  ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
	  return fluidViscosityScalar;
	}
	
	inline static ScalarProperty* createMinSeparationDistPb(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* minSeparationDistScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDist", caller);
	  return minSeparationDistScalar;
	}

	inline static ScalarProperty* createMaxSeparationDistRatioPb(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* maxSeparationDistRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistRatio", caller);
	  return maxSeparationDistRatioScalar;
	}
	inline static ScalarProperty* createFluidDensityPb(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* fluidDensityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidDensity", caller);
	  return fluidDensityScalar;
	}
}

namespace LIGGGHTS {
namespace ContactModels {

  template<>
  class CohesionModel<COHESION_ST_PF_LUBRICATION> : public CohesionModelBase {
  public:
	CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
		CohesionModelBase(lmp, hsetup, c),
		surfaceTension(0.0), contactAngle(NULL), fluidViscosity(0.0), minSeparationDist(0.),maxSeparationDistRatio(0.),fluidDensity(0.0)
	{
		
	}

	void registerSettings(Settings& settings) 
	{
		settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
		settings.registerOnOff("pressure",pressure_,true);
		settings.registerOnOff("capillary",capillary_,false);
	}

	inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

	void connectToProperties(PropertyRegistry & registry)
	{
		registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
		registry.registerProperty("contactAngle", &MODEL_PARAMS::createContactAngle);
		registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityPb);
		registry.registerProperty("minSeparationDist", &MODEL_PARAMS::createMinSeparationDistPb);
		registry.registerProperty("maxSeparationDistRatio", &MODEL_PARAMS::createMaxSeparationDistRatioPb);
		registry.registerProperty("fluidDensity", &MODEL_PARAMS::createFluidDensityPb);
		registry.connect("surfaceTension",surfaceTension ,"cohesion_model st/pf/lubrication");
		registry.connect("contactAngle", contactAngle,"cohesion_model st/pf/lubrication");
		registry.connect("fluidViscosity", fluidViscosity,"cohesion_model st/pf/lubrication");
		registry.connect("minSeparationDist", minSeparationDist,"cohesion_model st/pf/lubrication");
		registry.connect("maxSeparationDistRatio", maxSeparationDistRatio,"cohesion_model st/pf/lubrication");
		registry.connect("fluidDensity", fluidDensity,"cohesion_model st/pf/lubrication");
		// error checks on coarsegraining
		if(force->cg_active())
			error->cg(FLERR,"cohesion model st/pf/lubrication");

		neighbor->register_contact_dist_factor(maxSeparationDistRatio); 
		if(maxSeparationDistRatio < 1.0)
			error->one(FLERR,"\n\ncohesion model st/pf/lubrication requires maxSeparationDistanceRatio >= 1.0. Please increase this value.\n");
	}

	void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
	  if (!sidata.is_wall) { //r is the distance between the sphere's centers
	      const double ri = sidata.radi;
	      const double rj = sidata.radj;
	      const double zi = sidata.zi;
	      const double zj = sidata.zj;
	      const double deltan = sidata.deltan;
	      const double delx = sidata.delta[0];
	      const double delz = sidata.delta[2];
	      const double dxz = sqrt(delx*delx+delz*delz);
	      const double Ri = 1/(2*dxz)*sqrt((-dxz+ri-rj)*(-dxz-ri+rj)*(-dxz+ri+rj)*(dxz+ri+rj));
	      const double rhoi = sidata.densityi;
	      const double rhoj = sidata.densityj;
	  
	      double rb = 0, rp = 0, zp = 0 ,zb = 0;
	      if (rhoi > 0.1 && rhoj> 0.1) {        // make sure SI unit is used
		      rp = rhoi > 10 ? ri : rj;
		      rb = rhoi > 10 ? rj : ri;
		      zp = rhoi > 10 ? zi : zj;
		      zb = rhoi > 10 ? zi : zj;
	      } else {
		      rp = rhoi >= 1 ? ri : rj;
		      rb = rhoi >= 1 ? rj : ri;	
		      zp = rhoi >= 1 ? zi : zj;
		      zb = rhoi >= 1 ? zi : zj;
	      }
	      const double contactEffAngle = rp == rhoi ? contactAngle[sidata.itype] : contactAngle[sidata.itype];
	      //const double sin_alpha = Ri/rp;
	      /*const double hp = sqrt(rp * rp - Ri * Ri);
	      const double hb = sqrt(rb * rb - Ri * Ri);
	  
	      const double sp = 2*M_PI*rp*hp;
	      const double sb = 2*M_PI*rb*hb;

	      const double contactEffAngle = rp == rhoi ? contactAngle[sidata.itype] : contactAngle[sidata.itype];

	      const double wa = surfaceTension*(sb-sp*cos(contactEffAngle));

	      const double F_ad = -wa/deltan;
	      if(tangentialReduce_) sidata.Fn += F_ad; */
	      double Fca = 0,Fp = 0;
	      const double sinalpha = Ri/rb;
	      const double cosalpha = (dxz*dxz+rp*rp-rb*rb)/(2*dxz*rp);
	  
	      if (capillary_) {
		      const double sintheta_alpha = sin(contactEffAngle)*cosalpha-cos(contactEffAngle)*sinalpha;
		      Fca=-2*M_PI*surfaceTension*rp*sinalpha*sintheta_alpha;
	      }

	      if (pressure_) {
		      const double POb = sinalpha*rb;
		      const double ObM = abs(delz);
		      const double ObN = POb/dxz*ObM;
		      double H = 0;
		      if (zp - zb > 0) H = rb-ObN;
		      if (zp - zb > 0) H = rb+ObN;
		      Fp = M_PI*rp*rp*sinalpha*sinalpha*(2*surfaceTension/rb-fluidDensity*9.81*H);
	      }

	      const double fx = (Fca + Fp) * sidata.en[0];
	      const double fy = (Fca + Fp) * sidata.en[1];
	      const double fz = (Fca + Fp) * sidata.en[2];
	      if (tangentialReduce_) sidata.Fn += Fca + Fp;

	      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

		  i_forces.delta_F[0] += fx;
		  i_forces.delta_F[1] += fy;
		  i_forces.delta_F[2] += fz;

		  j_forces.delta_F[0] -= fx;
		  j_forces.delta_F[1] -= fy;
		  j_forces.delta_F[2] -= fz;
	    }
	}

	inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
	void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
	void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

	void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
	{
	   if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;

	   scdata.has_force_update = true;

	   if (!scdata.is_wall) {

		  double **v = atom->v;
		  const int i = scdata.i;
		  const int j = scdata.j;
		  const double rsq = scdata.rsq;
		  const double r = sqrt(rsq);
		  const double rinv =  1.0/r;
		  const double radsum = scdata.radsum;
		  const double radi = scdata.radi;
		  const double radj = scdata.radj;
		  const double rEff = radi*radj / radsum;
		  double d = r - radsum;
		  d = d > minSeparationDist ? d : minSeparationDist;
			  
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
		  const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
		  
		  const double F_lubrication = -6*M_PI*fluidViscosity*vn*rEff*rEff/d;
			  
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

	  if(scdata.is_wall) {
		
		double d = scdata.nonConDeltan;                             // deltan is the distance to the wall if scdata.wall = true
		d = d > minSeparationDist ? d : minSeparationDist;
		const double rinv =  1.0/scdata.nonConr;
		const double enx = scdata.delta[0] * rinv;
		const double eny = scdata.delta[1] * rinv;
		const double enz = scdata.delta[2] * rinv;
		const double rEff =  scdata.radi;
						
		const double vr1 = scdata.v_i[0]  - scdata.v_j[0];
		const double vr2 = scdata.v_i[1]  - scdata.v_j[1];
		const double vr3 = scdata.v_i[2]  - scdata.v_j[2];

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
	double surfaceTension;
	double * contactAngle;
	bool tangentialReduce_,capillary_,pressure_;
	double fluidViscosity,minSeparationDist, maxSeparationDistRatio,fluidDensity;
  };
}
}

#endif // COHESION_MODEL_ST_PF_LUBRICATION_H_
#endif
