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

	This cohesion model accounts for the particle-bubble adhesion force and force.
	Contributor: Linhan Ge
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_ST_PF,st/pf,9)
#else

#ifndef COHESION_MODEL_ST_PF_H_
#define COHESION_MODEL_ST_PF_H_

// #include <iostream>
#include "contact_models.h"
#include "cohesion_model_base.h"
#include "global_properties.h"
#include <cmath>
#include <math.h>
#include "neighbor.h"

namespace MODEL_PARAMS
{
	inline static MatrixProperty* createContactAngleHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  MatrixProperty* contactAngleMatrix = MODEL_PARAMS::createPerTypePairProperty(registry, "contactAngle", caller);
	  return contactAngleMatrix;
	}

    inline static MatrixProperty* createSurfaceTensionHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  MatrixProperty* surfaceTensionMatrix = MODEL_PARAMS::createPerTypePairProperty(registry, "surfaceTension", caller);
	  return surfaceTensionMatrix;
	}
}

namespace LIGGGHTS {
namespace ContactModels {

  template<>
  class CohesionModel<COHESION_ST_PF> : public CohesionModelBase {
  public:
	CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
		CohesionModelBase(lmp, hsetup, c),
		surfaceTension(NULL), 
		contactAngle(NULL), 
		liquidDensity(0.),
		tangentialReduce_(false),
		capillary_(true),
		pressure_(true)
	{
		
	}

	void registerSettings(Settings& settings) 
	{
		settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
		settings.registerOnOff("pressure",pressure_,true);
		settings.registerOnOff("capillary",capillary_,true);
	}

	inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

	void connectToProperties(PropertyRegistry & registry)
	{
		registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTensionHS);
		registry.registerProperty("contactAngle", &MODEL_PARAMS::createContactAngleHS);
		registry.registerProperty("liquidDensity", &MODEL_PARAMS::createLiquidDensity);
		registry.connect("surfaceTension",surfaceTension ,"cohesion_model st/pf");
		registry.connect("contactAngle", contactAngle,"cohesion_model st/pf");
		registry.connect("liquidDensity", liquidDensity,"cohesion_model st/pf");
		
		// error checks on coarsegraining
		if(force->cg_active())
			error->cg(FLERR,"cohesion model st/pf");
	}

	void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
	   
		//r is the distance between the sphere's centers
		const double ri = sidata.radi;
		const double rj = sidata.radj;
        const int itype = sidata.itype;
        const int jtype = sidata.jtype;
		const double zi = sidata.zi;
		const double zj = sidata.zj;

		const double rhoi = sidata.densityi;
		const double rhoj = sidata.densityj;
        const double deltan = sidata.deltan;
        
		double rb = 0, rp = 0, d = 0, Ri = 0, deltaz = 0, zp = 0, zb = 0;
        if ( rhoi/rhoj > 500 || rhoi/rhoj < 500 )
		{  
            rp = rhoi/rhoj > 500 ? ri : rj;
		    rb = rhoi/rhoj > 500 ? rj : ri;
		    zp = rhoi/rhoj > 500 ? zi : zj;
		    zb = rhoi/rhoj > 500 ? zi : zj;
            d = rb+rp-deltan;
            Ri = sqrt((-d+rp-rb)*(-d-rp+rb)*(-d+rp+rb)*(d+rp+rb))/(2*d);
			deltaz = abs(zi-zj);
        }
        else 
        {
            if(sidata.contact_flags) *sidata.contact_flags &= ~CONTACT_COHESION_MODEL;
            return;
        }

        const double theta = contactAngle[itype][jtype];
        const double sigma = surfaceTension[itype][jtype];
		  
		double Fca = 0,Fp = 0;
		const double sinalpha = Ri/rp;
		const double cosalpha = (d*d+rp*rp-rb*rb)/(2*d*rp);
	  
		if (capillary_) 
		{
		    const double sintheta_alpha = sin(theta)*cosalpha-cos(theta)*sinalpha;
			Fca = -2*M_PI*sigma*rp*sinalpha*sintheta_alpha;
			sidata.capillary = Fca;
            // std::cout << "capillary force " << Fca <<"  " << sinalpha << "  " << sintheta_alpha << std::endl; 
	    }

		if (pressure_) 
		{
		    const double n = rb - (rp-rb+d)*(rp+rb-d)/(2*d);
			const double m = deltaz/d*n;
            double H = 0;         
			if (zp > zb) H = rb - m;
			else  H = rb+m;
			Fp = M_PI*rp*rp*sinalpha*sinalpha*(liquidDensity*9.81*H-2*sigma/rb);
            // std::cout << "pressure force " << Fp <<"rb " << rb << std::endl; 
		}
         
        const double fx = (Fca + Fp) * sidata.en[0];
		const double fy = (Fca + Fp) * sidata.en[1];
		const double fz = (Fca + Fp) * sidata.en[2];

		if (tangentialReduce_) sidata.Fn += Fca + Fp;
          
		if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
        
        if (!sidata.is_wall) 
		{
            i_forces.delta_F[0] += fx;
		    i_forces.delta_F[1] += fy;
		    i_forces.delta_F[2] += fz;

		    j_forces.delta_F[0] -= fx;
		    j_forces.delta_F[1] -= fy;
		    j_forces.delta_F[2] -= fz;
		}
		else
		{
            i_forces.delta_F[0] = 0;
		    i_forces.delta_F[1] = 0;
		    i_forces.delta_F[2] = 0;

		    j_forces.delta_F[0] = 0;
		    j_forces.delta_F[1] = 0;
		    j_forces.delta_F[2] = 0;
		}
	}

	inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
	void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
	void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

	void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
	{
	    if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
	}

  private:
	double **surfaceTension;
	double **contactAngle;
	double liquidDensity;
	bool tangentialReduce_,capillary_,pressure_;
  };
}
}

#endif // COHESION_MODEL_ST_PF_LUBRICATION_H_
#endif
