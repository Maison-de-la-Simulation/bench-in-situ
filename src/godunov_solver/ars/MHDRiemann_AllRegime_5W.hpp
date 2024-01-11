#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace riemann
{

class MHD_AllRegime_5W
{
public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;

    MHD_AllRegime_5W(const Params& params, const EquationOfState& eos, const UniformGrid& grid);
    //MHD_AllRegime_5W(const MHD_AllRegime_5W& x) = default;
    //MHD_AllRegime_5W(MHD_AllRegime_5W&& x) = default;
    //~MHD_AllRegime_5W() = default;
    //MHD_AllRegime_5W& operator=(const MHD_AllRegime_5W& x) = default;
    //MHD_AllRegime_5W& operator=(MHD_AllRegime_5W&& x) = default;

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
MHD_AllRegime_5W::MHD_AllRegime_5W(const Params& params, const EquationOfState& eos, const UniformGrid& grid)
    : m_eos {eos}
    , m_K   {params.hydro.K}
    , m_g   {params.hydro.g}
    , m_dl  {grid.m_dl}
    , m_all_regime_correction {params.run.all_regime_correction}

{
}

KOKKOS_FORCEINLINE_FUNCTION
typename MHD_AllRegime_5W::ConsState
MHD_AllRegime_5W::operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const
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

    Vector<three_d, Real> ustar, pstar, cL, cR, clpcrm1, gdl, M;

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

      Real c = m_K*fmax(c_left, c_right);

      c_aL=c;
      c_aR=c;
      c_bR=c;
      c_bL=c;
    }
    cL(IX)   = c_aL; cR(IX)   = c_aR;
    cL(IY)   = c_aL; cR(IY)   = c_aR;
    cL(IZ)   = c_aL; cR(IZ)   = c_aR;
    cL(idir) = c_bL; cR(idir) = c_bR;

    clpcrm1(IX) = 1.0/(cL(IX)+cR(IX));
    clpcrm1(IY) = 1.0/(cL(IY)+cR(IY));
    clpcrm1(IZ) = 1.0/(cL(IZ)+cR(IZ));

    ustar(IX) = clpcrm1(IX)*(cL(IX)*q_L.v(IX) + cR(IX)*q_R.v(IX) + pL(IX) - pR(IX) + M(IX) ) ;
    ustar(IY) = clpcrm1(IY)*(cL(IY)*q_L.v(IY) + cR(IY)*q_R.v(IY) + pL(IY) - pR(IY) + M(IY) ) ;
    ustar(IZ) = clpcrm1(IZ)*(cL(IZ)*q_L.v(IZ) + cR(IZ)*q_R.v(IZ) + pL(IZ) - pR(IZ) + M(IZ) ) ;

    const Real Ma_norm = std::fabs(ustar(idir))/std::fmin(c0_L, c0_R);

    const Real theta_norm = m_all_regime_correction ? std::fmin(Ma_norm, 1.0) : 1.0;

    pstar(IX) = clpcrm1(IX)*(cR(IX)*pL(IX) + cL(IX)*pR(IX) + theta_norm*cL(IX)*cR(IX)*(q_L.v(IX) - q_R.v(IX)) ) ;
    pstar(IY) = clpcrm1(IY)*(cR(IY)*pL(IY) + cL(IY)*pR(IY) + theta_norm*cL(IY)*cR(IY)*(q_L.v(IY) - q_R.v(IY)) ) ;
    pstar(IZ) = clpcrm1(IZ)*(cR(IZ)*pL(IZ) + cL(IZ)*pR(IZ) + theta_norm*cL(IZ)*cR(IZ)*(q_L.v(IZ) - q_R.v(IZ)) ) ;

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
