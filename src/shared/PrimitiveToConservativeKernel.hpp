#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

class PrimitiveToConservativeKernel : public BaseKernel
{
    using Super = BaseKernel;

    using Euler = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using ConsState = typename Euler::ConsState;
    using PrimState = typename Euler::PrimState;

    using IntVector = IntVectorNd< three_d >;

public:
    PrimitiveToConservativeKernel( const Params& params, const UniformGrid& grid,
                                   const Array3d& q, const Array3d& u )
        : m_grid( grid )
        , m_u( u )
        , m_q( q )
        , m_eos( params.thermo )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( Int j ) const
    {
        const PrimState q_j = Super::getPrim( m_q, j );
        const ConsState u_j = Euler::primitiveToConservative( q_j, m_eos );
        Super::set( m_u, j, u_j );
    }

    const UniformGrid m_grid;
    const Array3d m_u;
    const Array3d m_q;
    const EquationOfState m_eos;
};

} // namespace hydro
