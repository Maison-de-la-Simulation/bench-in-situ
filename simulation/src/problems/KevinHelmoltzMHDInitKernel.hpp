#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include <cmath>

namespace hydro { namespace problems
{

class KevinHelmoltzMHDInitKernel : public BaseKernel
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
    KevinHelmoltzMHDInitKernel(const Params& params,
                        Array u, const UniformGrid& grid)
        : m_params {params}
        , m_eos    {params.thermo}
        , m_u      {u}
        , m_grid    {grid}
        , x_mid    {constants::half}
        , y_mid    {constants::half}
    {
    }

    static void apply(const Params& params,
                      Array u, const UniformGrid& grid)
    {
        KevinHelmoltzMHDInitKernel kernel {params, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("KevinHelmoltzMHD initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;
        const RealVector OX {m_grid.getCellCenter(j)};


        PrimState q_j;
        q_j.d     = 1.0;
        q_j.p     = 1.0/(m_eos.computeAdiabaticIndex());
        q_j.v(IX) = 0.5*tanh(20.0*OX[IY]);
        //q_j.v(IY) = 0.1*sin(2*constants::pi*OX[IX])*exp(-OX[IY]*OX[IY]/0.01);//papier edouard
        q_j.v(IY) = 0.01*sin(2*constants::pi*OX[IX])*exp(-OX[IY]*OX[IY]/0.01);//papier reference

        q_j.v(IZ) = 0.0;
        q_j.B(IX) = 0.1*cos(constants::pi/3);
        q_j.B(IY) = 0.0;
        q_j.B(IZ) = 0.1*sin(constants::pi/3);

        q_j.X=0.0;
        const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};

        set(m_u, j, u_j);
    }

    Params m_params;
    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
    Real x_mid, y_mid;
};

}}
