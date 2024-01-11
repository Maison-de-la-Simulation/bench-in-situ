#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"

namespace hydro { namespace problems
{

class BlastInitKernel : public BaseKernel
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
    BlastInitKernel(const Params& params, Array u, const UniformGrid& grid)
        : m_eos {params.thermo}
        , m_u      {u}
        , m_grid    {grid}
    {
    }

    static void apply(const Params& params, Array u,
                      const UniformGrid& grid)
    {
        BlastInitKernel kernel {params, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("Blast initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;

        const RealVector OX {m_grid.getCellCenter(j)};

        const Real r = sqrt(OX[IX]*OX[IX]+ OX[IY]*OX[IY]);

        PrimState q_j;

        q_j.d=1.0;

        q_j.v(IX) = 0.0;
        q_j.v(IY) = 0.0;
        q_j.v(IZ) = 0.0;

        q_j.B(IX) = sqrt(4*constants::pi)*1.0/sqrt(2.0);// High plasma beta
        q_j.B(IY) = sqrt(4*constants::pi)*1.0/sqrt(2.0);// High plasma beta

        q_j.B(IZ) = 0.0;

        if (r <0.1)
        {
          q_j.p = 10.0;// High plasma beta
        }
        else {
          q_j.p=0.1; 
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
