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
#include <cmath>

namespace MODEL_PARAMS {

	inline static ScalarProperty* createMaxSeparationDistRatioMHS(PropertyRegistry & registry, const char * caller, bool sanity_checks) // MHS means my hooke stiffness
	{
	  ScalarProperty* maxSeparationDistRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistRatio", caller);
	  return maxSeparationDistRatioScalar;
	}   

	inline static ScalarProperty* createFluidViscosityMHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
	  return fluidViscosityScalar;
	}

	inline static ScalarProperty* createFluidDensityMHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* fluidDensityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidDensity", caller);
	  return fluidDensityScalar;
	}
		inline static ScalarProperty* createEta_eMHS(PropertyRegistry & registry, const char * caller, bool sanity_checks)
	{
	  ScalarProperty* eta_eScalar = MODEL_PARAMS::createScalarProperty(registry, "eta_e", caller);
	  return eta_eScalar;
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
	  gamma_t(NULL),
	  e_dry(NULL),
	  eta_e(NULL),
	  tangential_damping(false),
	  limitForce(false),
	  displayedSettings(false),
	  history_offset(0),
	  impact_velocity(0),
	  elastic_potential_offset_(0),
	  elasticpotflag_(false),
	  fix_dissipated_(NULL),
	  dissipatedflag_(false),
	  dissipation_history_offset_(0)
	{
	  history_offset = hsetup->add_history_value("contflag", "0");
	  impact_velocity = hsetup->add_history_value("impactVelocity", "1");
	}

	void registerSettings(Settings & settings)
	{
	  settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("wallOnly", tangential_damping, false);
      settings.registerOnOff("limitForce", limitForce);
	  settings.registerOnOff("computeElasticPotential", elasticpotflag_, false);
	  settings.registerOnOff("computeDissipatedEnergy", dissipatedflag_, false);
	}

	inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
	{
		if (elasticpotflag_)
		{
			elastic_potential_offset_ = cmb->get_history_offset("elastic_potential_normal");
			if (elastic_potential_offset_ == -1)
			{
				elastic_potential_offset_ = hsetup->add_history_value("elastic_potential_normal", "0");
				hsetup->add_history_value("elastic_force_normal_0", "1");
				hsetup->add_history_value("elastic_force_normal_1", "1");
				hsetup->add_history_value("elastic_force_normal_2", "1");
				hsetup->add_history_value("elastic_torque_normal_i_0", "0");
				hsetup->add_history_value("elastic_torque_normal_i_1", "0");
				hsetup->add_history_value("elastic_torque_normal_i_2", "0");
				hsetup->add_history_value("elastic_torque_normal_j_0", "0");
				hsetup->add_history_value("elastic_torque_normal_j_1", "0");
				hsetup->add_history_value("elastic_torque_normal_j_2", "0");
				if (cmb->is_wall())
					hsetup->add_history_value("elastic_potential_wall", "0");
				cmb->add_history_offset("elastic_potential_normal", elastic_potential_offset_);
			}
		}
		if (dissipatedflag_)
		{
			if (cmb->is_wall())
			{
				fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_wall", "property/atom", "vector", 0, 0, "dissipated energy"));
				dissipation_history_offset_ = cmb->get_history_offset("dissipation_force");
				if (!dissipation_history_offset_)
					error->one(FLERR, "Internal error: Could not find dissipation history offset");
			}
			else
				fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy", "property/atom", "vector", 0, 0, "dissipated energy"));
			if (!fix_dissipated_)
				error->one(FLERR, "Surface model has not registered dissipated_energy fix");
		}

	}

	void connectToProperties(PropertyRegistry & registry) {
	  registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
	  registry.registerProperty("k_t", &MODEL_PARAMS::createKt);
	  registry.registerProperty("e_dry", &MODEL_PARAMS::createEdry);
	  registry.registerProperty("eta_e", &MODEL_PARAMS::createEta_eMHS);
	  registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityMHS);
	  registry.registerProperty("fluidDensity", &MODEL_PARAMS::createFluidDensityMHS);
	  registry.connect("k_n", k_n,"model myhooke/stiffness");
	  registry.connect("k_t", k_t,"model myhooke/stiffness");
	  registry.connect("e_dry", e_dry,"model myhooke/stiffness");
	  registry.connect("eta_e", eta_e,"model myhooke/stiffness");
	  registry.connect("fluidViscosity", fluidViscosity,"model myhooke/stiffness");
	  registry.connect("fluidDensity", fluidDensity,"model myhooke/stiffness");

      registry.registerProperty("gammat", &MODEL_PARAMS::createGammat);
	  registry.connect("gammat", gamma_t,"model myhooke/stiffness");

	  // error checks on coarsegraining
	  if(force->cg_active())
		error->cg(FLERR,"model myhooke/stiffness");

      neighbor->register_contact_dist_factor(1.05);

	  // enlarge contact distance flag in case of elastic energy computation
	  // to ensure that surfaceClose is called after a contact
	  if (elasticpotflag_)
		  //set neighbor contact_distance_factor here
		  neighbor->register_contact_dist_factor(1.01);
	}

	// effective exponent for stress-strain relationship
	
	inline double stressStrainExponent()
	{
	  return 1.;
	}

	void dissipateElasticPotential(SurfacesCloseData &scdata)
	{
		if (elasticpotflag_)
		{
			double * const elastic_energy = &scdata.contact_history[elastic_potential_offset_];
			if (scdata.is_wall)
			{
				// we need to calculate half an integration step which was left over to ensure no energy loss, but only for the elastic energy. The dissipation part is handled in fix_wall_gran_base.h.
				double delta[3];
				scdata.fix_mesh->triMesh()->get_global_vel(delta);
				vectorScalarMult3D(delta, update->dt);
				// -= because force is in opposite direction
				// no *dt as delta is v*dt of the contact position
				elastic_energy[0] -= (delta[0]*(elastic_energy[1]) +
									  delta[1]*(elastic_energy[2]) +
									  delta[2]*(elastic_energy[3]))*0.5
									 // from previous half step
									 + elastic_energy[10];
				elastic_energy[10] = 0.0;
			}
			elastic_energy[1] = 0.0;
			elastic_energy[2] = 0.0;
			elastic_energy[3] = 0.0;
			elastic_energy[4] = 0.0;
			elastic_energy[5] = 0.0;
			elastic_energy[6] = 0.0;
			elastic_energy[7] = 0.0;
			elastic_energy[8] = 0.0;
			elastic_energy[9] = 0.0;
		}
	}

	inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
	{
	  if (sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
	  
	  const int itype = sidata.itype;
	  const int jtype = sidata.jtype;
      const int i = sidata.i;
	  const int j = sidata.j;
	  const double rhoi = sidata.densityi;
	  const double rhoj = sidata.densityj;
	  const double radi = sidata.radi;
      const double r = sqrt(sidata.rsq);
      double meff=sidata.meff;

      double kn = k_n[itype][jtype];
	  double kt = k_t[itype][jtype];
      const double edry = e_dry[itype][jtype];
	  double gamman, gammat;
	  
	  if(!displayedSettings)
	  {
		displayedSettings = true;

		/*
		if(limitForce)
			if(0 == comm->me) fprintf(screen," NormalModel<MYHOOKE_STIFFNESS>: will limit normal force.\n");
		*/
	  }

	  gammat = gamma_t[itype][jtype];


	  if (!tangential_damping) gammat = 0.0;

	  // convert Kn and Kt from pressure units to force/distance^2
	  kn /= force->nktv2p;
	  kt /= force->nktv2p;
      
      // Eq.3.13 Izard, E., Bonometti, T., Lacaze, L., 2014. Modelling the dynamics of a sphere approaching and bouncing on a wall in a viscous fluid. Journal of Fluid Mechanics 747, 422-446.
      if (sidata.is_wall){
          double * const impactVelocity = &sidata.contact_history[impact_velocity]; 
	      const double  impactVn = impactVelocity[0];
          const double st = (rhoi +  0.5*fluidDensity)*impactVn*2*radi/fluidDensity;
	      const double stc = log(radi/eta_e);
          double ewet = 0;
          if (st < stc) {
              gamman = 2*sqrt(meff*kn);
          } 
          else {
              ewet = edry*(1-stc/st)*exp(-M_PI/2/sqrt(st-stc));     
              gamman = 2*log(ewet)*sqrt(meff*kn)/sqrt(log(ewet)*log(ewet) + M_PI*M_PI);
          }
      } else {
          if (wallOnly) gamman = 2*log(edry)*sqrt(meff*kn)/sqrt(log(edry)*log(edry) + M_PI*M_PI);
      }
      
      const double Fn_damping = -gamman*sidata.vn;    
	  const double Fn_contact = kn*sidata.deltan;
      double Fn = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force

      if(limitForce && (Fn<0.0) )
	  {
		  Fn = 0.0;
	  }

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
	  if (update_history)
	  {
		  // compute increment in elastic potential
		  if (elasticpotflag_)
		  {
			  double * const elastic_energy = &sidata.contact_history[elastic_potential_offset_];
			  // correct for wall influence
			  if (sidata.is_wall)
			  {
				  double delta[3];
				  sidata.fix_mesh->triMesh()->get_global_vel(delta);
				  vectorScalarMult3D(delta, update->dt);
				  // -= because force is in opposite direction
				  // no *dt as delta is v*dt of the contact position
					//printf("pela %e %e %e %e\n",  update->get_cur_time()-update->dt, deb, -sidata.radj, deb-sidata.radj);
				  elastic_energy[0] -= (delta[0]*elastic_energy[1] +
										delta[1]*elastic_energy[2] +
										delta[2]*elastic_energy[3])*0.5
									   // from previous half step
									   + elastic_energy[10];
				  elastic_energy[10] = -(delta[0]*Fn_contact*sidata.en[0] +
										 delta[1]*Fn_contact*sidata.en[1] +
										 delta[2]*Fn_contact*sidata.en[2])*0.5;
			  }
			  elastic_energy[1] = -Fn_contact*sidata.en[0];
			  elastic_energy[2] = -Fn_contact*sidata.en[1];
			  elastic_energy[3] = -Fn_contact*sidata.en[2];
			  elastic_energy[4] = 0.0;
			  elastic_energy[5] = 0.0;
			  elastic_energy[6] = 0.0;
			  elastic_energy[7] = 0.0;
			  elastic_energy[8] = 0.0;
			  elastic_energy[9] = 0.0;
		  }
		  // compute increment in dissipated energy
		  if (dissipatedflag_)
		  {
			  double * const * const dissipated = fix_dissipated_->array_atom;
			  double * const dissipated_i = dissipated[sidata.i];
			  double * const dissipated_j = dissipated[sidata.j];
			  const double F_diss = -Fn_damping;
			  dissipated_i[1] += sidata.en[0]*F_diss;
			  dissipated_i[2] += sidata.en[1]*F_diss;
			  dissipated_i[3] += sidata.en[2]*F_diss;
			  if (sidata.j < atom->nlocal && !sidata.is_wall)
			  {
				  dissipated_j[1] -= sidata.en[0]*F_diss;
				  dissipated_j[2] -= sidata.en[1]*F_diss;
				  dissipated_j[3] -= sidata.en[2]*F_diss;
			  }
			  else if (sidata.is_wall)
			  {
				  double * const diss_force = &sidata.contact_history[dissipation_history_offset_];
				  diss_force[0] -= sidata.en[0]*F_diss;
				  diss_force[1] -= sidata.en[1]*F_diss;
				  diss_force[2] -= sidata.en[2]*F_diss;
			  }
		  }
		  #ifdef NONSPHERICAL_ACTIVE_FLAG
		  if ((dissipatedflag_ || elasticpotflag_) && sidata.is_non_spherical)
			  error->one(FLERR, "Dissipation and elastic potential do not compute torque influence for nonspherical particles");
		  #endif
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
	  }
	}

	void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&)
	{
		if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
		
		double * const impactVelocity = &scdata.contact_history[impact_velocity]; 
		const int i = scdata.i;
		const int j = scdata.j;
		const double r = sqrt(scdata.rsq);
        const double radi = scdata.radi;
	    const double radj = scdata.is_wall ? radi : scdata.radj;
        			
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
		impactVelocity[0] = vn;
	}

	void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
	void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
	double ** k_n;
	double ** k_t;
	double ** gamma_t;
    double ** e_dry;

	bool tangential_damping;
	bool limitForce,wallOnly;
	bool displayedSettings;
	int elastic_potential_offset_;
	bool elasticpotflag_;
	FixPropertyAtom *fix_dissipated_;
	bool dissipatedflag_;
	int dissipation_history_offset_;
	int history_offset,impact_velocity;
	double eta_e,fluidViscosity,fluidDensity;
  };
}
}
#endif // NORMAL_MODEL_MYHOOKE_STIFFNESS_H_
#endif