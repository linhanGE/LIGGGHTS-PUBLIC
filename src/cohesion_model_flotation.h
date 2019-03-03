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

	Linhan Ge (University of Newcastle)

	Copyright 2012-     DCS Computing GmbH, Linz
	Copyright 2009-2012 JKU Linz

	This cohesion model accounts for the particle-bubble adhesion force and force.
	Contributor: Linhan Ge
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_FLOTATION,flotation,9)
#else

#ifndef COHESION_MODEL_FLOTATION
#define COHESION_MODEL_FLOTATION

// #include <iostream>
#include "contact_models.h"
#include "cohesion_model_base.h"
#include "global_properties.h"
#include <cmath>
#include <math.h>
#include "neighbor.h"

namespace MODEL_PARAMS
{
	inline static MatrixProperty* createContactAngleF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  MatrixProperty* contactAngleMatrix = MODEL_PARAMS::createPerTypePairProperty(registry, "contactAngle", caller);
	  return contactAngleMatrix;
	}

    inline static MatrixProperty* createSurfaceTensionF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  MatrixProperty* surfaceTensionMatrix = MODEL_PARAMS::createPerTypePairProperty(registry, "surfaceTension", caller);
	  return surfaceTensionMatrix;
	}

	inline static ScalarProperty* createXDLVOCutOffF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
		ScalarProperty* XDLVOCutOff_Scalar = MODEL_PARAMS::createScalarProperty(registry, "XDLVOCutOff", caller);
		return XDLVOCutOff_Scalar;
	}

	inline static ScalarProperty* createAF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
		ScalarProperty* A_Scalar = MODEL_PARAMS::createScalarProperty(registry, "A", caller);
		return A_Scalar;
	}

	inline static ScalarProperty* createKF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
		ScalarProperty* K_Scalar = MODEL_PARAMS::createScalarProperty(registry, "K", caller);
		return K_Scalar;
	}

    inline static ScalarProperty* createLamdaF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
		ScalarProperty* lamda_Scalar = MODEL_PARAMS::createScalarProperty(registry, "lamda", caller);
		return lamda_Scalar;
	}

	inline static MatrixProperty* createMinSeparationDistF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  MatrixProperty* minSeparationDistMatrix = MODEL_PARAMS::createPerTypePairProperty(registry, "minSeparationDist", caller);
	  return minSeparationDistMatrix;
	}
    
	inline static ScalarProperty* createMaxSeparationDistRatioF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* maxSeparationDistRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistRatio", caller);
	  return maxSeparationDistRatioScalar;
	}

	inline static ScalarProperty* createMuF(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* muScalar = MODEL_PARAMS::createScalarProperty(registry, "mu", caller);
	  return muScalar;
	}
}

namespace LIGGGHTS {
namespace ContactModels {

  template<>
  class CohesionModel<COHESION_FLOTATION> : public CohesionModelBase {
  public:
	CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
		CohesionModelBase(lmp, hsetup, c),
		surfaceTension(NULL), 
		contactAngle(NULL), 
		XDLVOCutOff(0.),
		A(0.),
		K(0.),
		lamda(0.),
		tangentialReduce_(false),
		capillary_(true),
		lubrication_(true),
		vdw_(true),
		hydrophobic_(true)
	{
		
	}

	void registerSettings(Settings& settings) 
	{
		settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
		settings.registerOnOff("capillary",capillary_,true);
		settings.registerOnOff("lubrication",lubrication_,true);
		settings.registerOnOff("vdw",vdw_,true);
		settings.registerOnOff("hydrophobic",hydrophobic_,true);
	}

	inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

	void connectToProperties(PropertyRegistry & registry)
	{
		registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTensionF);
		registry.registerProperty("contactAngle", &MODEL_PARAMS::createContactAngleF);
		registry.registerProperty("XDLVOCutOff", &MODEL_PARAMS::createXDLVOCutOffF);
		registry.registerProperty("A", &MODEL_PARAMS::createAF);
		registry.registerProperty("K", &MODEL_PARAMS::createKF);
        registry.registerProperty("lamda", &MODEL_PARAMS::createLamdaF);
		registry.registerProperty("mu", &MODEL_PARAMS::createMuF);
		registry.registerProperty("minSeparationDist", &MODEL_PARAMS::createMinSeparationDistF);
		registry.registerProperty("maxSeparationDistRatio", &MODEL_PARAMS::createMaxSeparationDistRatioF);
                
		registry.connect("surfaceTension",surfaceTension ,"cohesion_model flotation");
		registry.connect("contactAngle", contactAngle,"cohesion_model flotation");
        registry.connect("XDLVOCutOff", XDLVOCutOff, "cohesion_model flotation");
		registry.connect("A", A, "cohesion_model flotation");
		registry.connect("K",K ,"cohesion_model flotation");
        registry.connect("lamda",lamda ,"cohesion_model flotation");
		registry.connect("mu", mu,"cohesion_model flotation");
		registry.connect("minSeparationDist", minSeparationDist,"cohesion_model flotation");
		registry.connect("maxSeparationDistRatio", maxSeparationDistRatio,"cohesion_model flotation");
		
		neighbor->register_contact_dist_factor(1.1);
	}

	void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
		//r is the distance between the sphere's centers
		const double ri = sidata.radi;
		const double rj = sidata.radj;
        const int itype = sidata.itype;
        const int jtype = sidata.jtype;

		const double rhoi = sidata.densityi;
		const double rhoj = sidata.densityj;
        const double d = sqrt(sidata.rsq);
		
		const double theta = contactAngle[itype][jtype];
        const double sigma = surfaceTension[itype][jtype];

		double rb = 0, rp = 0, Ri = 0, Fca=0, sinalpha=0, cosalpha=0, sintheta_alpha=0;

        if ( rhoi/rhoj > 500 || rhoi/rhoj < 500 )
		{  
	        if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

            rp = rhoi/rhoj > 500 ? ri : rj;
		    rb = rhoi/rhoj > 500 ? rj : ri;
            Ri = sqrt((-d+rp-rb)*(-d-rp+rb)*(-d+rp+rb)*(d+rp+rb))/(2*d);
			sinalpha = Ri/rp;
		    cosalpha = (d*d+rp*rp-rb*rb)/(2*d*rp);
	        sintheta_alpha = sin(theta)*cosalpha-cos(theta)*sinalpha;
        }
        else 
        {
            if(sidata.contact_flags) *sidata.contact_flags &= ~CONTACT_COHESION_MODEL;
                return;
        }
		
		if (capillary_) 
		{
     		Fca = -2*M_PI*sigma*Ri*sintheta_alpha;
			sidata.capillary = Fca;
            // std::cout << "capillary force " << Fca <<"  " << sinalpha << "  " << sintheta_alpha << std::endl; 
	    }
        
        const double fx = Fca * sidata.en[0];
		const double fy = Fca * sidata.en[1];
		const double fz = Fca * sidata.en[2];
        
		if (tangentialReduce_) sidata.Fn += Fca;
       
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
    	const int itype = scdata.itype;
		const int jtype = scdata.jtype;
		const double rhoi = scdata.densityi;
		const double rhoj = scdata.densityj;
		const double rsq = scdata.rsq;
		const double r = sqrt(rsq);
		const double rinv =  1.0/r;
		const double radsum = scdata.radsum;
		const double radi = scdata.radi;
		const double radj = scdata.radj;
		const double rEff = radi*radj / radsum;
	    
	    if (rhoi/rhoj > 500 || rhoi/rhoj < 500)
			if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
		else 
			return;
		
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
		
        double F_lubrication=0;

		if (d <= maxSeparationDistRatio*radsum && lubrication_)
		    F_lubrication = -6*M_PI*mu*vn*rEff*rEff/d;
        
		double Fvdw = 0, Fh = 0;
		double H = r - radsum;
		H = std::max(XDLVOCutOff,H);

        // hydrophobic force
		if ( hydrophobic_) Fh = -rEff*K*exp(-H/lamda);
        
		// van-der walls force
		if (vdw_) Fvdw = A*rEff/(6*H*H);

		const double fx = (F_lubrication + Fvdw + Fh) * enx;       //en represent the normal direction vector, en[0] is the x coordinate
		const double fy = (F_lubrication + Fvdw + Fh) * eny;				 
		const double fz = (F_lubrication + Fvdw + Fh) * enz;		
		
		scdata.has_force_update = true;
		
        if (!scdata.is_wall)
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

  private:
	double **surfaceTension;
	double **contactAngle;
    double XDLVOCutOff,A,K,lamda;	
	double mu, maxSeparationDistRatio;
	double ** minSeparationDist;
	bool tangentialReduce_,capillary_,lubrication_,vdw_,hydrophobic_;
  };
}
}

#endif // COHESION_MODEL_ST_PF_LUBRICATION_H_
#endif
