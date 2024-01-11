#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace riemann
{

class MHD3W_optimized
{
public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;

    MHD3W_optimized(const Params& params, const EquationOfState& eos, const UniformGrid& grid);

    KOKKOS_FORCEINLINE_FUNCTION ConsState
    operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const;

  private:
    const EquationOfState m_eos;
    const Real m_dz;
    const bool m_all_regime_correction;
    const Real m_gz;
    const Real m_inv_gamma_m1;

};
inline
MHD3W_optimized::MHD3W_optimized(const Params& params, const EquationOfState& eos, const UniformGrid& grid)
    : m_eos {eos}
    , m_dz  {grid.m_dl[IZ]}
    , m_all_regime_correction {params.run.all_regime_correction}
    , m_gz  {params.hydro.g[IZ]}
    , m_inv_gamma_m1{eos.computeInvGammam1()}
{
}

KOKKOS_FORCEINLINE_FUNCTION
typename MHD3W_optimized::ConsState
MHD3W_optimized::operator()(const PrimState &q_L, const PrimState &q_R, const int& side,const int& idir, bool st_powell) const
{
    using namespace constants;

    Vector<three_d, Real> pL = -q_L.B( idir )*q_L.B;
    pL(idir) += q_L.p + MHD::computeMagneticEnergy(q_L);

    Vector<three_d, Real> pR = -q_R.B( idir )*q_R.B;
    pR(idir) += q_R.p + MHD::computeMagneticEnergy(q_R);

    const Real M = 0.5*(q_L.d + q_R.d)*m_dz*m_gz;

    const Real a  = 1.1 * std::fmax(q_L.d*MHD::computeFastMagnetoAcousticSpeed(q_L, m_eos, idir),
                                    q_R.d*MHD::computeFastMagnetoAcousticSpeed(q_R, m_eos, idir));

    Vector<three_d, Real> ustar = 0.5*(q_R.v+q_L.v)-0.5*(1.0/a)*(pR-pL);
    ustar(IZ) += 0.5*(1.0/a)*M;

    Real theta = std::fabs(ustar(idir))/std::fmin(MHD::computeSpeedOfSound(q_L, m_eos)
                                           , MHD::computeSpeedOfSound(q_R, m_eos));

    theta = m_all_regime_correction ? std::fmin(theta, 1.0) : 1.0;

    //if (theta<0.02){theta=0.02;} // This was only for Grand Challenge

    const Vector<three_d, Real> pstar = 0.5*(pR+pL)-0.5*a*theta*(q_R.v-q_L.v);

    Real B_next;
    PrimState q;
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

     // if (st_powell)
     // {
     //   if (side==-1)
     //   {
     //     B_next = q_R.B(idir);
     //   }
     //   else
     //   {
     //     B_next = q_L.B(idir);
     //   }
     // }

    ConsState flux;

    flux.d = ustar(idir)*q.d;
    flux.m = ustar(idir)*q.d*q.v +pstar ;
    flux.e = ustar(idir)*(q.p*m_inv_gamma_m1 + 0.5*q.d*dot(q.v,q.v) +0.5*dot(q.B,q.B))+ dot(pstar,ustar);
    flux.B = ustar(idir)*q.B - B_next*ustar;
    flux.dX = ustar(idir)*q.d*q.X;

    if (idir==IZ)
    {
      flux.m ( IZ) += -side*M*0.5;
      flux.e += - 0.5*side*M*ustar(IZ);
    }

    return flux;
}

}}
