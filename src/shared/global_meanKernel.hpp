#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "global_mean.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

class global_meanKernel : public BaseKernel
{
    using Super = BaseKernel;

    using EquationOfState = MHDSystem::EquationOfState;
    using ConsState = MHDSystem::ConsState;
    using PrimState = MHDSystem::PrimState;
    using VC = MHDSystem::VarCons;
    using VP = MHDSystem::VarPrim;

public:
    global_meanKernel( const Params&, const UniformGrid& grid, Array3d q )
        : m_grid( grid )
        , m_q( q )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void init( global_mean& dst ) const
    {
        using constants::zero;
        dst.m_emag = zero;
        dst.m_ekin = zero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( Int i, Int j, Int k, global_mean& sum ) const
    {
        const PrimState q_j = Super::getPrim( m_q, m_grid.coordToIndex({i,j,k}) );
        // get primitive variables in current cell

        sum.m_emag += MHDSystem::computeMagneticEnergy(q_j);
        sum.m_ekin += MHDSystem::computeKineticEnergy(q_j);
    }

    KOKKOS_INLINE_FUNCTION
    void join( global_mean& dst, const global_mean& src ) const
    {
        // sum reduce
        dst.m_emag += src.m_emag;
        dst.m_ekin += src.m_ekin;
    }

    KOKKOS_INLINE_FUNCTION
    void join( volatile global_mean& dst, const volatile global_mean& src ) const
    {
        // sum reduce
        dst.m_emag += src.m_emag;
        dst.m_ekin += src.m_ekin;
    }

    UniformGrid m_grid;
    ConstArray3d m_q;
};

} // namespace hydro
