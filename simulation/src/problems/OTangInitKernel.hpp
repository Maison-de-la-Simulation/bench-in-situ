#pragma once

#include "OTangParams.hpp"
#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"

#include <cmath>

namespace hydro { namespace problems
{

class OTangInitKernel : public BaseKernel
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
    OTangInitKernel(const Params& params, const OTangParams& pParams,
                     Array u, const UniformGrid& grid)
        : m_params {params}
        , m_problem_params {pParams}
        , m_eos {params.thermo}
        , m_u {u}
        , m_grid {grid}
    {
    }

    static void apply(const Params& params, const OTangParams& pParams,
                      Array u, const UniformGrid& grid)
    {
        OTangInitKernel kernel {params, pParams, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("OTang initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;

        const RealVector OX {m_grid.getCellCenter(j)};

        PrimState q_j;
        const Real B0       = 1.0/(sqrt(4.0*constants::pi));
        q_j.d = 25.0/(36.0*constants::pi);
        q_j.p = 5.0 /(12.0*constants::pi);
        q_j.v( IX ) = -std::sin(2*constants::pi * OX[IY]);
        q_j.v( IY ) =  std::sin(2*constants::pi * OX[IX]);
        q_j.v( IZ ) = 0.0;
        q_j.B( IX ) = -B0*std::sin(2.0*constants::pi * OX[IY]);
        q_j.B( IY ) =  B0*std::sin(4.0*constants::pi * OX[IX]);
        q_j.B( IZ ) = 0.0;

        q_j.X=0.0;
        const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
        set(m_u, j, u_j);
    }

    Params m_params;
    OTangParams m_problem_params;
    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
};

}}
