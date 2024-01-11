#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace riemann
{

class MHD_AllRegime_3W
{
public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;

    MHD_AllRegime_3W(const Params& params, const EquationOfState& eos, const UniformGrid& grid);
    //MHD_AllRegime_3W(const MHD_AllRegime_3W& x) = default;
    //MHD_AllRegime_3W(MHD_AllRegime_3W&& x) = default;
    //~MHD_AllRegime_3W() = default;
    //MHD_AllRegime_3W& operator=(const MHD_AllRegime_3W& x) = default;
    //MHD_AllRegime_3W& operator=(MHD_AllRegime_3W&& x) = default;

    KOKKOS_FORCEINLINE_FUNCTION ConsState
    operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const;

  private:
    const EquationOfState m_eos;
    const Real m_K;
    const RealVector3d m_g;
    const RealVector3d m_dl;
    const bool m_all_regime_correction;
    const bool m_entropic_correction;


};
inline
MHD_AllRegime_3W::MHD_AllRegime_3W(const Params& params, const EquationOfState& eos, const UniformGrid& grid)
    : m_eos {eos}
    , m_K   {params.hydro.K}
    , m_g   {params.hydro.g}
    , m_dl  {grid.m_dl}
    , m_all_regime_correction {params.run.all_regime_correction}
    , m_entropic_correction {params.run.entropic_correction}


{
}

KOKKOS_FORCEINLINE_FUNCTION
typename MHD_AllRegime_3W::ConsState
MHD_AllRegime_3W::operator()(const PrimState &q_L, const PrimState &q_R, const int& side,const int& idir, bool st_powell) const
{
    using namespace constants;

    const Real rho_L = q_L.d;
    const Real p_L   = q_L.p;
          Real emag  = MHD::computeMagneticEnergy(q_L);
    const Real c_L   = rho_L*MHD::computeFastMagnetoAcousticSpeed(q_L, m_eos, idir);
    const Real c0_L   = MHD::computeSpeedOfSound(q_L, m_eos);

    Vector<three_d, Real> pL = -q_L.B( idir )*q_L.B;
    pL(idir) += p_L + emag;

    const Real rho_R = q_R.d;
    const Real p_R   = q_R.p;
               emag  = MHD::computeMagneticEnergy(q_R);
    const Real c_R   = rho_R*MHD::computeFastMagnetoAcousticSpeed(q_R, m_eos, idir);
    const Real c0_R   = MHD::computeSpeedOfSound(q_R, m_eos);

    Vector<three_d, Real> pR = -q_R.B( idir )*q_R.B;
    pR(idir) += p_R + emag;

    const Vector<three_d, Real> gdl(m_g[IX]*m_dl[IX], m_g[IY]*m_dl[IY], m_g[IZ]*m_dl[IZ]);

    const Vector<three_d, Real> M = 0.5*(q_L.d + q_R.d)*gdl;

    const Real a    = m_K * std::fmax(c_R, c_L);
    const Real a_m1 = 1.0/a;

    const Vector<three_d, Real> ustar = 0.5*(q_R.v+q_L.v)-0.5*a_m1*(pR-pL-M);

    const Real Ma_norm = std::fabs(ustar(idir))/std::fmin(c0_L, c0_R);
    Real theta_norm;
    if (m_entropic_correction)
    {
      const Real du = q_R.v(idir)-q_L.v(idir)+10.0e-15;
      const Real dp = (q_R.p-q_L.p)/a;

      const Real Am=du/(du-dp);
      const Real Ap=du/(du+dp);

      Real theta_m;
      Real theta_p;

      if (Am>0.0)
      {
        theta_m=dp/du;
      }
      else
      {
        theta_m=0;
      }

      if (Ap>0.0)
      {
        theta_p=-dp/du;
      }
      else
      {
        theta_p=0;
      }
      theta_norm = m_all_regime_correction ? fmax(fmin(fmax(theta_m,theta_p),1.),Ma_norm) : 1.0;

    //  if (theta_norm>Ma_norm){
    //    std::cout<<"yo "<<idir<<"  "<<theta_norm<<" "<<Ma_norm<<std::endl;
    //  }
    }
    else
    {
      theta_norm = m_all_regime_correction ? std::fmin(Ma_norm, 1.0) : 1.0;
    }

    const Vector<three_d, Real> pstar = 0.5*(pR+pL)-0.5*a*theta_norm*(q_R.v-q_L.v);

    PrimState q;
    Real B_next;

    if ( ustar(idir) > zero )
     {
       q = q_L;
       B_next = q_R.B(idir);
     }
     else
     {
       q = q_R;
       B_next = q_L.B(idir);
     }

     //Powell source term

     if (st_powell)
     {
       if (side==-1)
       {
         B_next = q_R.B(idir);
       }
       else
       {
         B_next = q_L.B(idir);
       }
     }

     const ConsState u = MHD::primitiveToConservative(q, m_eos);

     ConsState flux;

      flux.d = ustar(idir)*u.d;
      flux.m = ustar(idir)*u.m +pstar ;
      flux.m ( idir) += -side*M(idir)*0.5;
      flux.e = ustar(idir)*u.e + dot(pstar,ustar);
      flux.e += - 0.5*side*M(idir)*ustar(idir);
      flux.B = ustar(idir)*u.B - B_next*ustar;
      flux.dX = ustar(idir)*u.dX;

    return flux;
}

}}
