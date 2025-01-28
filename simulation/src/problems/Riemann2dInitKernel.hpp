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

class Riemann2dInitKernel : public BaseKernel
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
    Riemann2dInitKernel(const Params& params,
                        Array u, const UniformGrid& grid)
        : m_params {params}
        , m_eos    {params.thermo}
        , m_u      {u}
        , m_grid    {grid}
        , x_mid    {0.8}
        , y_mid    {0.8}
    {
    }

    static void apply(const Params& params,
                      Array u, const UniformGrid& grid)
    {
        Riemann2dInitKernel kernel {params, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("Riemann2d initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;
        const RealVector OX {m_grid.getCellCenter(j)};

        PrimState q_j;
        if (OX[IX]>x_mid)
        {
            if (OX[IY]>y_mid)
            {
                q_j.d = 1.5;
                q_j.p = 1.5;
                q_j.v( IX ) = 0.0;
                q_j.v( IY ) = 0.0;
                q_j.v( IZ ) = zero;
                q_j.B( IX ) = zero;
                q_j.B( IY ) = zero;
                q_j.B( IZ ) = zero;

            }
            else
            {
              // quarter 0
              q_j.d = 0.5323;
              q_j.p = 0.3;
              q_j.v( IX ) = 0.0;
              q_j.v( IY ) = 1.206;
              q_j.v( IZ ) = zero;
              q_j.B( IX ) = zero;
              q_j.B( IY ) = zero;
              q_j.B( IZ ) = zero;

            }
        }
        else
        {
            if (OX[IY]>y_mid)
            {
                // quarter 2
                q_j.d = 0.5323;
                q_j.p = 0.3;
                q_j.v( IX ) = 1.206;
                q_j.v( IY ) = 0.0;
                q_j.v( IZ ) = zero;
                q_j.B( IX ) = zero;
                q_j.B( IY ) = zero;
                q_j.B( IZ ) = zero;

            }
            else
            {
              // quarter 1
              q_j.d = 0.138;
              q_j.p = 0.029;
              q_j.v( IX ) = 1.206;
              q_j.v( IY ) = 1.206;
              q_j.v( IZ ) = zero;
              q_j.B( IX ) = zero;
              q_j.B( IY ) = zero;
              q_j.B( IZ ) = zero;


            }
        }

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
