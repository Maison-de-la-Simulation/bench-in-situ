#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include "slice_average.hpp"
#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

class slice_averageKernel : public BaseKernel
{
    using Super = BaseKernel;

    using Euler = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using ConsState = typename Euler::ConsState;
    using PrimState = typename Euler::PrimState;
    using VC = typename Euler::VarCons;
    using VP = typename Euler::VarPrim;

public:
    slice_averageKernel( const Params& params, const UniformGrid& grid, Array3d q, Int iz )
        : m_grid( grid )
        , m_eos( params.thermo )
        , m_q( q )
        , m_iz( iz )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void init( slice_average& dst ) const noexcept
    {
        dst.m_Bx2 = 0.0;
        dst.m_By2 = 0.0;
        dst.m_Bz2 = 0.0;
        dst.m_u = 0.0;
        dst.m_v = 0.0;
        dst.m_logmu = 0.0;
        dst.m_logtheta = 0.0;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const Int i, const Int j, slice_average& slice ) const noexcept
    {
        using namespace constants;
        const PrimState q_j = Super::getPrim( m_q, m_grid.coordToIndex({i, j, m_iz}) );

        slice.m_Bx2 += q_j.B(IX) * q_j.B(IX);
        slice.m_By2 += q_j.B(IY) * q_j.B(IY);
        slice.m_Bz2 += q_j.B(IZ) * q_j.B(IZ);
        slice.m_u += q_j.v(IX);
        slice.m_v += q_j.v(IY);
        slice.m_logmu += std::log(MHDSystem::computeMeanMolecularWeight(q_j, m_eos));
        slice.m_logtheta += MHDSystem::computeLogPotentialTemperature(q_j, m_eos);
    }

    KOKKOS_INLINE_FUNCTION
    void join( slice_average& dst, const slice_average& src ) const noexcept
    {
        dst.m_Bx2 += src.m_Bx2;
        dst.m_By2 += src.m_By2;
        dst.m_Bz2 += src.m_Bz2;
        dst.m_u += src.m_u;
        dst.m_v += src.m_v;
        dst.m_logmu += src.m_logmu;
        dst.m_logtheta += src.m_logtheta;
    }

    KOKKOS_INLINE_FUNCTION
    void join( volatile slice_average& dst, const volatile slice_average& src ) const noexcept
    {
        dst.m_Bx2 += src.m_Bx2;
        dst.m_By2 += src.m_By2;
        dst.m_Bz2 += src.m_Bz2;
        dst.m_u += src.m_u;
        dst.m_v += src.m_v;
        dst.m_logmu += src.m_logmu;
        dst.m_logtheta += src.m_logtheta;
    }

    UniformGrid m_grid;
    EquationOfState m_eos;
    ConstArray3d m_q;
    Int m_iz;
};

} // namespace hydro
