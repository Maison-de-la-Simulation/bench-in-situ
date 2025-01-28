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

class ConservativeToPrimitiveKernel : public BaseKernel
{
    using Super = BaseKernel;

    using Euler = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using ConsState = typename Euler::ConsState;
    using PrimState = typename Euler::PrimState;

public:
    ConservativeToPrimitiveKernel( const Params& params, const UniformGrid& grid,
                                   const Array3d& u, const Array3d& q )
        : m_grid( grid )
        , m_q( q )
        , m_u( u )
        , m_eos( params.thermo )
    {
    }

    //ConservativeToPrimitiveKernel( ConservativeToPrimitiveKernel const& x ) = default;

    //ConservativeToPrimitiveKernel( ConservativeToPrimitiveKernel&& x ) = default;

    ~ConservativeToPrimitiveKernel() = default;

    //ConservativeToPrimitiveKernel& operator=( ConservativeToPrimitiveKernel const& x ) = default;

    //ConservativeToPrimitiveKernel& operator=( ConservativeToPrimitiveKernel&& x ) = default;

    KOKKOS_INLINE_FUNCTION
    void operator()( Int j ) const
    {
        const ConsState u_j = Super::getCons( m_u, j );
        const PrimState q_j = Euler::conservativeToPrimitive( u_j, m_eos );
        Super::set( m_q, j, q_j );
    }

    const UniformGrid m_grid;
    const Array3d m_q;
    const ConstArray3d m_u;
    const EquationOfState m_eos;
};

} // namespace hydro
