#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"

namespace hydro { namespace problems
{

class ImplodeInitKernel : public BaseKernel
{
    using Super           = BaseKernel;

    using Euler           = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using ConsState       = typename Euler::ConsState;
    using PrimState       = typename Euler::PrimState;
    using RealVector      = RealVector3d;
    using VC              = typename Euler::VarCons;
    using VP              = typename Euler::VarPrim;
    static constexpr int nbvar = Euler::nbvar;

    using Array = Array3d;

public:
    ImplodeInitKernel(const Params& params, Array u, const UniformGrid& grid)
        : m_eos {params.thermo}
        , m_u      {u}
        , m_grid    {grid}
    {
    }

    static void apply(const Params& params, Array u,
                      const UniformGrid& grid)
    {
        ImplodeInitKernel kernel {params, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("Implode initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;

        const RealVector OX {m_grid.getCellCenter(j)};
        Real tmp {0.0};
        for (Int idim=0; idim<three_d; ++idim)
        {
            tmp += OX[idim];
        }

        PrimState q_j;
        q_j.d = tmp < half ? one : eigth;
        q_j.p = tmp < half ? one : tenth+four*hundredth;
        for (int idim=0; idim<three_d; ++idim)
        {
            q_j.v(idim) = 0.0;
            q_j.B(idim) = 0.0;
        }

        q_j.X=0.0;
        const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
        Super::set(m_u, j, u_j);
    }

    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
};

}}
