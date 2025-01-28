#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"

namespace hydro { namespace problems
{

class RotorInitKernel : public BaseKernel
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
    RotorInitKernel(const Params& params, Array u, const UniformGrid& grid)
        : m_eos {params.thermo}
        , m_u      {u}
        , m_grid    {grid}
    {
    }

    static void apply(const Params& params, Array u,
                      const UniformGrid& grid)
    {
        RotorInitKernel kernel {params, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("Rotor initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;

        const RealVector OX {m_grid.getCellCenter(j)};

        const Real x = OX[IX];
        const Real y = OX[IY];

        const Real r = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
        const Real r0 = 0.1;
        const Real r1 = 0.115;
        const Real U0 = 2.0;

        const Real f = (r1-r)/(r1-r0);

        PrimState q_j;

        q_j.p=1.0;
        q_j.B(IX) = 5.0/sqrt(4*constants::pi);
        q_j.B(IY) = 0.0;
        q_j.B(IZ) = 0.0;
        q_j.v(IZ) = 0.0;

        if(r<r0)
        {
          q_j.d = 10;
          q_j.v(IX) = -U0*(y-0.5)/r0;
          q_j.v(IY) =  U0*(x-0.5)/r0;
        }
        else if ((r<r1)&&(r>r0))
        {
          q_j.d = 1.0+9*f;
          q_j.v(IX) = -f*U0*(y-0.5)/r;
          q_j.v(IY) =  f*U0*(x-0.5)/r;
        }
        else
        {
          q_j.d = 1;
          q_j.v(IX) = 0.0;
          q_j.v(IY) = 0.0;
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
