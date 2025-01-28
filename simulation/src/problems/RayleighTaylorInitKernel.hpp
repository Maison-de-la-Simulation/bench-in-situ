#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "RayleighTaylorParams.hpp"
#include <cmath>

namespace hydro { namespace problems
{

class RayleighTaylorInitKernel : BaseKernel
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
    RayleighTaylorInitKernel(const Params& params, const RayleighTaylorParams& pParams,
                             Array u, const UniformGrid& grid)
        : m_params {params}
        , m_problem_params {pParams}
        , m_eos            {params.thermo}
        , m_u {u}
        , m_grid {grid}
    {
    }

    static void apply(const Params& params, const RayleighTaylorParams& pParams,
                      Array u, const UniformGrid& grid)
    {
        RayleighTaylorInitKernel kernel {params, pParams, u, grid};

        Int start_x = 0;
        Int start_z = 0;

        Int end_x = grid.m_nbCells[IX]+2*grid.m_ghostWidths[IX];
        Int end_z = grid.m_nbCells[IZ]+2*grid.m_ghostWidths[IZ];

        Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy ({start_x, start_z}, {end_x, end_z});
        Kokkos::parallel_for(policy, kernel);

    }

    KOKKOS_INLINE_FUNCTION
    Real phi(Int index) const
    {
        const RealVector OX{m_grid.getCellCenter(index)};
        Real phi{constants::zero};
        for (Int idim = 0; idim < two_d; ++idim)
        {
            phi -= m_params.hydro.g[idim] * OX[idim];
        }
        return phi;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int ix, Int iz) const
    {

        using namespace constants;

        for (Int iy=m_grid.m_ghostWidths[IY];
             iy<m_grid.m_nbCells[IY]+m_grid.m_ghostWidths[IY]; ++iy)
        {


            const IntVector coord {{ix, iy, iz}};
            const Int j {m_grid.coordToIndex(coord)};
            const IntVector coord_k {{ix, iy-1, iz}};
            const Int k {m_grid.coordToIndex(coord_k)};
            const ConsState u_k {getCons(m_u, k)};
            const PrimState q_k {Euler::conservativeToPrimitive(u_k, m_eos)};
            const RealVector OX {m_grid.getCellCenter(coord)};
            const Real Lx = m_params.mesh.up[IX] - m_params.mesh.low[IX];
            const Real Ly = m_params.mesh.up[IY] - m_params.mesh.low[IY];

            PrimState q_j;
            q_j.d = OX[IY] <= 0.0 ? m_problem_params.density_bottom : m_problem_params.density_up;
            q_j.p = iy == m_grid.m_ghostWidths[IY] ? m_problem_params.pressure_bottom : q_k.p + half*(q_k.d + q_j.d)*(phi(k) - phi(j));
            q_j.v( IX ) = 0.0;
            q_j.v( IZ ) = 0.0;
            q_j.B( IX ) = 0.0;
            q_j.B( IY ) = 0.0;
            q_j.B( IZ ) = 0.0;

            q_j.v( IY ) = m_problem_params.perturbation * (one+std::cos(pi_2*OX[IX]/Lx))*(one+std::cos(pi_2*OX[IY]/Ly))/four;

            q_j.X=0.0;
            const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
            set(m_u, j, u_j);
        }
    }

    Params m_params;
    RayleighTaylorParams m_problem_params;
    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
};

}}
