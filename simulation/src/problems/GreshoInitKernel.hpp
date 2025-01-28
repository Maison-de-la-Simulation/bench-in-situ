#pragma once

#include "GreshoParams.hpp"
#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace problems
{

class GreshoInitKernel : public BaseKernel
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
    GreshoInitKernel(const Params& params, const GreshoParams& pParams,
                     Array u, const UniformGrid& grid)
        : m_params {params}
        , m_problem_params {pParams}
        , m_eos {params.thermo}
        , m_u {u}
        , m_grid {grid}
    {
    }

    static void apply(const Params& params, const GreshoParams& pParams,
                      Array u, const UniformGrid& grid)
    {
        GreshoInitKernel kernel {params, pParams, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("Gresho initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;

        const RealVector OX {m_grid.getCellCenter(j)};
        const RealVector CX {{OX[IX]-m_problem_params.OC[IX], OX[IY]-m_problem_params.OC[IY]}};
        const Real r {std::sqrt(CX[IX]*CX[IX] + CX[IY]*CX[IY])};
        const Real theta {std::atan2(CX[IY], CX[IX])};

        Real vtheta {zero};
        Real p {zero};
        if (r <= two*tenth)
        {
            vtheta = five * r;
            p      = m_problem_params.p0 + (five*five*half) * r * r;
        }
        else if (r <= four*tenth)
        {
            vtheta = two - five * r;
            p      = m_problem_params.p0 + (five*five*half) * r * r + four * (one - five * r + std::log(five*r));
        }
        else
        {
            vtheta = zero;
            p      = m_problem_params.p0 - two + std::log(four*four);
        }

        PrimState q_j;
        q_j.d = m_problem_params.rho;
        q_j.p = p;
        q_j.v( IX ) = - vtheta * std::sin(theta);
        q_j.v( IY ) = + vtheta * std::cos(theta);
        q_j.v( IZ ) = + vtheta * std::cos(theta);
        q_j.B( IX ) = zero;
        q_j.B( IY ) = zero;
        q_j.B( IZ ) = zero;

        q_j.X=0.0;
        const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
        set(m_u, j, u_j);
    }

    Params m_params;
    GreshoParams m_problem_params;
    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
};

}}
