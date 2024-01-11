#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace riemann
{

class AllRegime_exp
{
public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;

    AllRegime_exp(const Params& params, const EquationOfState& eos, const UniformGrid& grid);
    //AllRegime_exp(const AllRegime_exp& x) = default;
    //AllRegime_exp(AllRegime_exp&& x) = default;
    //~AllRegime_exp() = default;
    //AllRegime_exp& operator=(const AllRegime_exp& x) = default;
    //AllRegime_exp& operator=(AllRegime_exp&& x) = default;

    KOKKOS_FORCEINLINE_FUNCTION ConsState
    operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const;//, const Real& dt) const

  private:
    const EquationOfState m_eos;
    const Real m_K;
    const RealVector3d m_g;
    const RealVector3d m_dl;
    const bool m_all_regime_correction;
    const Int m_all_regime_correction_type;


};
inline
AllRegime_exp::AllRegime_exp(const Params& params, const EquationOfState& eos, const UniformGrid& grid)
    : m_eos {eos}
    , m_K   {params.hydro.K}
    , m_g   {params.hydro.g}
    , m_dl  {grid.m_dl}
    , m_all_regime_correction {params.run.all_regime_correction}
    , m_all_regime_correction_type {params.run.all_regime_correction_type}
{
}

KOKKOS_FORCEINLINE_FUNCTION
typename AllRegime_exp::ConsState
AllRegime_exp::operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell ) const//, const Real& dt) const
{
    using namespace constants;

    const Real un_L = q_L.v( idir );
    const Real un_R = q_R.v( idir );

    // Find the largest eigenvalues in the normal direction to the interface
    const Real c_L = MHD::computeSpeedOfSound(q_L, m_eos);
    const Real c_R = MHD::computeSpeedOfSound(q_R, m_eos);

    const Real a = m_K * std::fmax(q_R.d * c_R, q_L.d * c_L);

    const Real M = 0.5*(q_L.d + q_R.d)*m_g[idir]*m_dl[idir];
    const Real ustar = half * (un_R+un_L) - half * (q_R.p-q_L.p - M) / a;

    const Real Ma = std::fabs(ustar)/(0.5*(c_L+c_R));

    const Int strictly_all_regime = 0;
    const Int chalons_2014_entropy_satisfying = 1;
    const Int semi_discrete_acoustic_entropy_satisfying = 2;
    const Int discrete_acoustic_entropy_satisfying = 3;

    Real theta = 1.0;

    Real y;
    Real Delu;

    Real pistar;

    if (m_all_regime_correction)
    {
      if (m_all_regime_correction_type == strictly_all_regime)
      {
        theta = Ma;
      }
      else if (m_all_regime_correction_type==chalons_2014_entropy_satisfying)
      {
        const Real tau_L  =  1.0/q_L.d;
        const Real tau_Ls = tau_L + (1.0/a) * ( ustar-un_L );
        const Real p_inter_L = q_L.p*pow(tau_L / tau_Ls, m_eos.computeAdiabaticIndex());

        const Real tau_R  =  1.0/q_R.d;
        const Real tau_Rs = tau_R - (1.0/a) * ( ustar-un_R  );
        const Real p_inter_R = q_R.p*pow(tau_R / tau_Rs, m_eos.computeAdiabaticIndex());

        const Real pistar_base = half * (q_R.p + q_L.p) - half * a * (un_R - un_L);

        const Real grad_p_L = (pistar_base - p_inter_L);
        const Real grad_p_R = (pistar_base - p_inter_R);

        const Real gradu = (un_R-un_L);
        const Real alpha_L = fabs(grad_p_L/gradu)*2/a;
        const Real alpha_R = fabs(grad_p_R/gradu)*2/a;

        theta = 1.0 - fmin(alpha_L,alpha_R);
        theta = fmin(1.0, theta);
        theta = fmax(0.0, theta);

      }
      else if (m_all_regime_correction_type==semi_discrete_acoustic_entropy_satisfying)
      {

        const Real du = un_R-un_L+10.0e-15;
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

        theta = fmax(fmin(fmax(theta_m,theta_p),1.),Ma);

      }
      else if (m_all_regime_correction_type==discrete_acoustic_entropy_satisfying)
      {
        Delu = un_R-un_L;
        const Real DelA = (un_R - q_R.p/a) - (un_L - q_L.p/a);
        const Real DelB = (un_R + q_R.p/a) - (un_L + q_L.p/a);

        const Real l_L = 0.0;//dt/(q_L.d*m_dl[idir]);
        const Real l_R = 0.0;//dt/(q_R.d*m_dl[idir]);

        const Real m2 = -( (l_L*a-0.5)*DelA*DelA + (l_R*a-0.5)*DelB*DelB );
        const Real q  = Delu - l_L*a*DelA - l_R*a*DelB;
        const Real Delta = q*q + 2*a*m2*(l_R+l_L);

        const Real sqDel = sqrt(Delta);

        const Real yp = ( -(Delu - l_L*a*DelA - l_R*a*DelB) + sqDel ) /(a*(l_L+l_R));
        const Real ym = ( -(Delu - l_L*a*DelA - l_R*a*DelB) - sqDel ) /(a*(l_L+l_R));


        if (Delu > 0)
        {
          y=fmin(yp, Delu);
        }
        else
        {
          y=fmax(ym,Delu);
        }

        const Real du = un_R-un_L;
        const Real dp = (q_R.p-q_L.p)/a;

        const Real Am = du/(du-dp);
        const Real Ap = du/(du+dp);

        Real thetaSD;
        if (Am <= 0.0) {} else {thetaSD =  dp/du;}
        if (Ap <= 0.0) {} else {thetaSD = -dp/du;}

        Real ySD=(1.0-thetaSD)*du;

      //  if((l_L*a-0.5>0)||(l_R*a-0.5>0)){std::cout<<"chelou"<<std::endl;}

      }
    }

    pistar = half * (q_R.p + q_L.p) - half * theta* a * (un_R - un_L);

    if( m_all_regime_correction_type == discrete_acoustic_entropy_satisfying)
    {
      pistar = half * (q_R.p + q_L.p)-(a/2)*(Delu-y);
    }

    PrimState q;
    if ( ustar > zero )
    {
        q = q_L;
    }
    else
    {
        q = q_R;
    }
    ConsState u = MHD::primitiveToConservative(q,m_eos);

    ConsState flux;

    flux.d = ustar*u.d;
    flux.m = ustar*u.m;
    flux.e = ustar*u.e;
    flux.dX= ustar*u.dX;
    flux.B = 0.0*q.B;

    flux.m(idir) +=  pistar - side*M*0.5;
    flux.e +=  pistar*ustar - side*M*ustar*0.5;

    if (st_powell){;}//Do nothing, powell is only for MHD
    return flux;
}

}}
