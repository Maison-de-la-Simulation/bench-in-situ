#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "TimeStep.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

class TimeStepKernel : public BaseKernel
{
    using Super = BaseKernel;

    using Euler = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using ConsState = typename Euler::ConsState;
    using PrimState = typename Euler::PrimState;
    using RealVector = RealVector3d;
    using VC = typename Euler::VarCons;
    using VP = typename Euler::VarPrim;
    static constexpr int nbvar = Euler::nbvar;

public:
    TimeStepKernel( const Params& params, const UniformGrid& grid, Array3d q )
        : m_grid( grid )
        , m_eos( params.thermo )
        , m_q( q )
        , m_invdl {}
        , m_magnetic_resistivity_enabled( params.hydro.magnetic_resistivity_enabled )
        , m_resitivity_coefficient( params.hydro.resistivity_coefficient )

    {
        for ( int idim = 0; idim < three_d; ++idim )
        {
            m_invdl[ idim ] = constants::one / grid.m_dl[ idim ];
        }
    }

    KOKKOS_INLINE_FUNCTION
    void init( TimeStep& dst ) const
    {
        using constants::zero;
        dst.hydro = zero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( Int j, TimeStep& invDt ) const
    {
        using namespace constants;

        const PrimState q_j = Super::getPrim( m_q, j );
        // get primitive variables in current cell


        Real invDtLoc = zero;

        const Real B_norm22 = dot(q_j.B, q_j.B);
        const Real c    = m_eos.computeSpeedOfSound(q_j.d, q_j.p);
        const Real c02  = c*c;
        const Real ca2  = B_norm22/q_j.d;
        const Real cap2x = q_j.B(IX)*q_j.B(IX)/q_j.d;
        const Real cap2y = q_j.B(IY)*q_j.B(IY)/q_j.d;
        const Real cap2z = q_j.B(IZ)*q_j.B(IZ)/q_j.d;

        const Real c_jx = sqrt(0.5*(c02+ca2)+0.5*sqrt((c02+ca2)*(c02+ca2)-4.*c02*cap2x));
        const Real c_jy = sqrt(0.5*(c02+ca2)+0.5*sqrt((c02+ca2)*(c02+ca2)-4.*c02*cap2y));
        const Real c_jz = sqrt(0.5*(c02+ca2)+0.5*sqrt((c02+ca2)*(c02+ca2)-4.*c02*cap2z));


        const Real u_jx = q_j.v( IX );
        const Real u_jy = q_j.v( IY );
        const Real u_jz = q_j.v( IZ );

        invDtLoc += m_invdl[ IX ] * ( c_jx + std::fabs( u_jx ) );
        invDtLoc += m_invdl[ IY ] * ( c_jy + std::fabs( u_jy ) );
        invDtLoc += m_invdl[ IZ ] * ( c_jz + std::fabs( u_jz ) );

        if (m_magnetic_resistivity_enabled)

        { Real invDt_diffusion = 2*m_resitivity_coefficient * ( m_invdl[ IX ] * m_invdl [ IX ] + m_invdl[ IY ] * m_invdl [ IY ] + m_invdl[ IZ ] * m_invdl [ IZ ] );

          invDtLoc = std::fmax(invDt_diffusion, invDtLoc);

        }

        if ( invDt.hydro < invDtLoc )
        {
            invDt.hydro = invDtLoc;
        }


    }

    KOKKOS_INLINE_FUNCTION
    void join(  TimeStep& dst, const  TimeStep& src ) const
    {
        // max reduce
        if ( dst.hydro < src.hydro )
        {
            dst.hydro = src.hydro;
        }
    }

    KOKKOS_INLINE_FUNCTION
    void join( volatile TimeStep& dst, const volatile TimeStep& src ) const
    {
        // max reduce
        if ( dst.hydro < src.hydro )
        {
            dst.hydro = src.hydro;
        }
    }

    const UniformGrid m_grid;
    const EquationOfState m_eos;
    const ConstArray3d m_q;
    RealVector m_invdl;
    const bool m_magnetic_resistivity_enabled;
    const Real m_resitivity_coefficient;

};

} // namespace hydro
