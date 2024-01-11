#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace riemann
{

class AllRegime
{
public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;

    AllRegime(const Params& params, const EquationOfState& eos, const UniformGrid& grid);
    //AllRegime(const AllRegime& x) = default;
    //AllRegime(AllRegime&& x) = default;
    //~AllRegime() = default;
    //AllRegime& operator=(const AllRegime& x) = default;
    //AllRegime& operator=(AllRegime&& x) = default;

    KOKKOS_FORCEINLINE_FUNCTION ConsState
    operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const;

  private:
    const EquationOfState m_eos;
    const Real m_K;
    const RealVector3d m_g;
    const RealVector3d m_dl;
    const bool m_all_regime_correction;

};
inline
AllRegime::AllRegime(const Params& params, const EquationOfState& eos, const UniformGrid& grid)
    : m_eos {eos}
    , m_K   {params.hydro.K}
    , m_g   {params.hydro.g}
    , m_dl  {grid.m_dl}
    , m_all_regime_correction {params.run.all_regime_correction}

{
}

KOKKOS_FORCEINLINE_FUNCTION
typename AllRegime::ConsState
AllRegime::operator()(const PrimState &q_L, const PrimState &q_R, const int& side, const int& idir, bool st_powell) const
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
    const Real Ma = std::fmax(std::fabs(un_L)/c_L, std::fabs(un_R)/c_R);
    const Real theta = m_all_regime_correction ? std::fmin(Ma, one) : 1.0;
    const Real pistar = half * (q_R.p + q_L.p) - half * theta * a * (un_R - un_L);

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
