#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "RiemannParams.hpp"

namespace hydro { namespace problems
{

class RiemannInitKernel : public BaseKernel
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

    using Array           = Array3d;
    using IntVector       = IntVectorNd<three_d>;

public:
    RiemannInitKernel(const Params& params, const RiemannParams& pParams,
                      Array u, const UniformGrid& grid)
        : m_eos            {params.thermo}
        , m_u              {u}
        , m_grid           {grid}
        , m_problem_params {pParams}
    {
    }

    static void apply(const Params& params, const RiemannParams& pParams,
                      Array u, const UniformGrid& grid)
    {
        RiemannInitKernel kernel {params, pParams, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("Riemann initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        const RealVector OX {m_grid.getCellCenter(j)};
        Int IXd=0, IYd=0, IZd=0, dir = m_problem_params.direction;
        if ( dir == IX)
        {
          IXd = IX; IYd = IY; IZd = IZ;
        }
        else if (dir == IY)
        {
          IXd = IY; IYd = IZ; IZd = IX;
        }
        else if (dir == IZ)
        {
          IXd = IZ; IYd = IX; IZd = IY;
        }

        PrimState q_j;
        q_j.d = OX[dir] < m_problem_params.center ? m_problem_params.density_left : m_problem_params.density_right;
        q_j.p = OX[dir] < m_problem_params.center ? m_problem_params.pressure_left : m_problem_params.pressure_right;
        q_j.v( IXd ) = OX[dir] < m_problem_params.center ? m_problem_params.x_velocity_left : m_problem_params.x_velocity_right;
        q_j.v( IYd ) = OX[dir] < m_problem_params.center ? m_problem_params.y_velocity_left : m_problem_params.y_velocity_right;
        q_j.v( IZd ) = OX[dir] < m_problem_params.center ? m_problem_params.z_velocity_left : m_problem_params.z_velocity_right;

        q_j.B( IXd ) = OX[dir] < m_problem_params.center ? m_problem_params.x_Mag_Field_left : m_problem_params.x_Mag_Field_right;
        q_j.B( IYd ) = OX[dir] < m_problem_params.center ? m_problem_params.y_Mag_Field_left : m_problem_params.y_Mag_Field_right;
        q_j.B( IZd ) = OX[dir] < m_problem_params.center ? m_problem_params.z_Mag_Field_left : m_problem_params.z_Mag_Field_right;

        q_j.X=0.0;
        const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
        Super::set(m_u, j, u_j);
    }

    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
    RiemannParams m_problem_params;
};

}}
