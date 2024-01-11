#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace riemann
{

class MHD_ARWB_5W
{
public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;

    MHD_ARWB_5W(const Params& params, const EquationOfState& eos, const UniformGrid& grid);
    //MHD_ARWB_5W(const MHD_ARWB_5W& x) = default;
    //MHD_ARWB_5W(MHD_ARWB_5W&& x) = default;
    //~MHD_ARWB_5W() = default;
    //MHD_ARWB_5W& operator=(const MHD_ARWB_5W& x) = default;
    //MHD_ARWB_5W& operator=(MHD_ARWB_5W&& x) = default;

    KOKKOS_FORCEINLINE_FUNCTION ConsState
    operator()(const PrimState &q_L, const PrimState &q_R, const int& side,const int& idir, bool st_powell) const;

  private:
    const EquationOfState m_eos;
    const Real m_K;
    const RealVector3d m_g;
    const RealVector3d m_dl;
    const bool m_all_regime_correction;

};
inline
MHD_ARWB_5W::MHD_ARWB_5W(const Params& params, const EquationOfState& eos, const UniformGrid& grid)
    : m_eos {eos}
    , m_K   {params.hydro.K}
    , m_g   {params.hydro.g}
    , m_dl  {grid.m_dl}
    , m_all_regime_correction {params.run.all_regime_correction}

{
}

KOKKOS_FORCEINLINE_FUNCTION
typename MHD_ARWB_5W::ConsState
MHD_ARWB_5W::operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const
{
    using namespace constants;

    const Real rho_L       = q_L.d;
    const Real p_L         = q_L.p;
          Real emag        = MHD::computeMagneticEnergy(q_L);
          Real B_norm2     = q_L.B(idir)*q_L.B(idir);
          Real B_trans2    = 2*emag - B_norm2;
          Real BnormBtrans = sqrt(B_norm2*B_trans2);
          Real c0_L          = MHD::computeSpeedOfSound(q_L,m_eos);
          Real c_aL        = sqrt(rho_L*(B_norm2 + BnormBtrans))+1e-14;
          Real c_bL        = sqrt(rho_L*(rho_L*c0_L*c0_L + B_trans2 + BnormBtrans));

    Vector<three_d, Real> pL;
    pL        = -q_L.B( idir )*q_L.B;
    pL(idir) += p_L + emag;

    const Real rho_R       = q_R.d;
    const Real p_R         = q_R.p;
               emag        = MHD::computeMagneticEnergy(q_R);
               B_norm2     = q_R.B(idir)*q_R.B(idir);
               B_trans2    = 2*emag - B_norm2;
               BnormBtrans = sqrt(B_norm2*B_trans2);
          Real c0_R        = MHD::computeSpeedOfSound(q_R,m_eos);
          Real c_aR        = sqrt(rho_R*(B_norm2 + BnormBtrans))+1e-14;
          Real c_bR        = sqrt(rho_R*(rho_R*c0_R*c0_R + (B_trans2 + BnormBtrans)));

    Vector<three_d, Real> pR;
    pR        = -q_R.B( idir )*q_R.B;
    pR(idir) += p_R + emag;

    Vector<three_d, Real> ustar, pstar, c, denom, gdl, M;

    Real c_a, c_b;

    gdl(IX) = m_g[IX]*m_dl[IX];
    gdl(IY) = m_g[IY]*m_dl[IY];
    gdl(IZ) = m_g[IZ]*m_dl[IZ];

    M = 0.5*(q_L.d + q_R.d)*gdl;

    if((q_L.B(IX)*q_R.B(IX)<-1e-16)
    || (q_L.B(IY)*q_R.B(IY)<-1e-16)
    || (q_L.B(IZ)*q_R.B(IZ)<-1e-16))
    {
      const Real c_left  = rho_L*MHD::computeFastMagnetoAcousticSpeed(q_L, m_eos, idir);
      const Real c_right = rho_R*MHD::computeFastMagnetoAcousticSpeed(q_R, m_eos, idir);

      Real cK = m_K*fmax(c_left, c_right);

      c_a=cK;
      c_b=cK;
    }
    else
    {
      c_a = m_K*fmax(c_aL, c_aR);
      c_b = m_K*fmax(c_bL, c_bR);
    }
    c(IX)   = c_a;
    c(IY)   = c_a;
    c(IZ)   = c_a;
    c(idir) = c_b;

    denom(IX) = 1.0/(2*c(IX));
    denom(IY) = 1.0/(2*c(IY));
    denom(IZ) = 1.0/(2*c(IZ));

    ustar(IX) = denom(IX)*(c(IX)*q_L.v(IX) + c(IX)*q_R.v(IX) + pL(IX) - pR(IX) + M(IX) ) ;
    ustar(IY) = denom(IY)*(c(IY)*q_L.v(IY) + c(IY)*q_R.v(IY) + pL(IY) - pR(IY) + M(IY) ) ;
    ustar(IZ) = denom(IZ)*(c(IZ)*q_L.v(IZ) + c(IZ)*q_R.v(IZ) + pL(IZ) - pR(IZ) + M(IZ) ) ;

    const Real Ma_norm = std::fabs(ustar(idir))/std::fmin(c0_L, c0_R);

    const Real theta_norm = m_all_regime_correction ? std::fmin(Ma_norm, 1.0) : 1.0;

    pstar(IX) = denom(IX)*(c(IX)*pL(IX) + c(IX)*pR(IX) + theta_norm*c(IX)*c(IX)*(q_L.v(IX) - q_R.v(IX)) ) ;
    pstar(IY) = denom(IY)*(c(IY)*pL(IY) + c(IY)*pR(IY) + theta_norm*c(IY)*c(IY)*(q_L.v(IY) - q_R.v(IY)) ) ;
    pstar(IZ) = denom(IZ)*(c(IZ)*pL(IZ) + c(IZ)*pR(IZ) + theta_norm*c(IZ)*c(IZ)*(q_L.v(IZ) - q_R.v(IZ)) ) ;

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

     ConsState u = MHD::primitiveToConservative(q, m_eos);

     ConsState flux;

     flux.d = ustar(idir)*u.d;
     flux.m = ustar(idir)*u.m +pstar;
     flux.m ( idir) += -side*M(idir)*0.5;
     flux.e = ustar(idir)*u.e + dot(pstar,ustar);
     flux.e += - 0.5*side*M(idir)*ustar(idir);
     flux.B = ustar(idir)*u.B - B_next*ustar;
     flux.dX = ustar(idir)*u.dX;
    return flux;

}

}}
