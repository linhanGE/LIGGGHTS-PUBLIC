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

#ifdef NORMAL_MODEL
NORMAL_MODEL(MYHOOKE_STIFFNESS,myhooke/stiffness,7)
#else
#ifndef NORMAL_MODEL_MYHOOKE_STIFFNESS_H_
#define NORMAL_MODEL_MYHOOKE_STIFFNESS_H_
#include "contact_models.h"
#include "normal_model_base.h"

namespace MODEL_PARAMS {
	inline static ScalarProperty* createTcMHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
		ScalarProperty* tc_Scalar = MODEL_PARAMS::createScalarProperty(registry, "tc", caller);
		return tc_Scalar;
	}

    inline static ScalarProperty* createBubbleIDMHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
		ScalarProperty* bubbleID_Scalar = MODEL_PARAMS::createScalarProperty(registry, "bubbleID", caller);
		return bubbleID_Scalar;
	}
}

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<MYHOOKE_STIFFNESS> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
      NormalModelBase(lmp, hsetup, c),
      k_n(NULL),
      k_t(NULL),
      betaeff(NULL),
      tc(0.),
      bubbleID(1.),
      collision_offset_(0),
      collisionflag_(false),
      tangential_damping(false),
      fully_damping(false),
      elasticForceOff(false),
      limitForce(false),
      displayedSettings(false)
    {
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("fully_damping", fully_damping, false);
      settings.registerOnOff("elasticForceOff", elasticForceOff, false);
      settings.registerOnOff("limitForce", limitForce);
      settings.registerOnOff("computeCollision", collisionflag_, false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
        if (collisionflag_)
        {
            collision_offset_ = hsetup->add_history_value("collisionCount","0");
            hsetup->add_history_value("collisionIndicator", "0");
            hsetup->add_history_value("impactVelocity", "0");
            hsetup->add_history_value("collisionIDIndicator", "0");			
        }
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
      registry.registerProperty("k_t", &MODEL_PARAMS::createKt);
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff);
      registry.registerProperty("tc", &MODEL_PARAMS::createTcMHS);
      registry.registerProperty("bubbleID", &MODEL_PARAMS::createBubbleIDMHS);

      registry.connect("k_n", k_n,"model myhooke/stiffness");
      registry.connect("k_t", k_t,"model myhooke/stiffness");
      registry.connect("betaeff", betaeff,"model myhooke/stiffness");
      registry.connect("tc", tc, "model myhooke/stiffness");
      registry.connect("bubbleID", bubbleID, "model myhooke/stiffness");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model myhooke/stiffness");
      
      // enlarge contact distance flag in case of elastic energy computation
      // to ensure that surfaceClose is called after a contact
      if (collisionflag_)
          //set neighbor contact_distance_factor here
          neighbor->register_contact_dist_factor(1.1);
    }

	// effective exponent for stress-strain relationship
    inline double stressStrainExponent()
    {
      return 1.;
    } 
      
    void collisionCompute(SurfacesCloseData &scdata)
    {
        if (collisionflag_)
        {
            double * const collision_time = &scdata.contact_history[collision_offset_];
            
            collision_time[0] = 0.0;
            collision_time[1] = 0.0;
            collision_time[2] = 0.0;
            collision_time[3] = 0.0;
        }
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      
      if (sidata.contact_flags)
          *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
                
      const bool update_history = sidata.computeflag && sidata.shearupdate;

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const int i = sidata.i;
      const int j = sidata.j;

      double meff=sidata.meff;  // already consider the freeze particle in pair_gran_base.h Line 394

      double kn = k_n[itype][jtype];
      double kt = k_t[itype][jtype];
      double gamman = 0;
      
      if (fully_damping)
          gamman = 2*sqrt(meff*kn);
      else 
          gamman = -2*sqrt(meff*kn)*betaeff[itype][jtype]; // betaeff is negative, gamman should be positive

      const double gammat = tangential_damping ? gamman : 0.0;

      if(!displayedSettings)
      {
        displayedSettings = true;
      }

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;
            
      double Fn, Fn_damping, Fn_contact;

      if (elasticForceOff)
          Fn = 0;
      else
      {
          Fn_damping = -gamman*sidata.vn;    
          Fn_contact = kn*sidata.deltan;
          Fn = Fn_damping + Fn_contact;    
      }
      
      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          double torque_i[3] = {0.0, 0.0, 0.0}; //initialized here with zeros to avoid compiler warnings
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
      #endif

      // energy balance terms
      if (update_history && !elasticForceOff)
      {
          if (collisionflag_)
          {
              double * const collision_time = &sidata.contact_history[collision_offset_];
              // correct for wall influence
              if (MathExtraLiggghts::compDouble(collision_time[0], 0, 1e-16)) 
              {
                  collision_time[1] = 1;
                  collision_time[2] = sidata.vn;
                  if ( atom->tag[i] == int (bubbleID) || atom->tag[j] == int (bubbleID) )
                      collision_time[3] = 1;
                  else 
                      collision_time[3] = 0;
              }
              else
              {
                  collision_time[1] = 0;
                  collision_time[3] = 0;
              } 

              collision_time[0] += 1;
		      
			  if ( collision_time[0] >= tc)
                  sidata.fluidContactNormal = 0;
	          else sidata.fluidContactNormal = 1;
          }
      }

      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];
                }
        #endif
      } else {
        i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  double xcj[3], torque_j[3];
                  double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
                  vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
                  vectorCross3D(xcj, Fn_j, torque_j);

                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];

                  j_forces.delta_torque[0] += torque_j[0];
                  j_forces.delta_torque[1] += torque_j[1];
                  j_forces.delta_torque[2] += torque_j[2];
                }
        #endif
      }
    }

    void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&)
    {   
        if (scdata.contact_flags)
            *scdata.contact_flags |= CONTACT_NORMAL_MODEL;
        
        collisionCompute(scdata);
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** k_n;
    double ** k_t;
    double ** betaeff;
    double tc, bubbleID;
    int collision_offset_;
    bool collisionflag_;
    bool tangential_damping, fully_damping, elasticForceOff;
    bool limitForce;
    bool displayedSettings;
  };
}
}
#endif // NORMAL_MODEL_MYHOOKE_STIFFNESS_H_
#endif
