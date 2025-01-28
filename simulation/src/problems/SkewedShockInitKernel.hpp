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

class SkewedShockInitKernel : public BaseKernel
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
    SkewedShockInitKernel(const Params& params,
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
        SkewedShockInitKernel kernel {params, u, grid};
        Kokkos::RangePolicy<Int> range {0, grid.nbCells()};
        Kokkos::parallel_for("SkewedShock initialization kernel - RangePolicy", range, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(Int j) const
    {
        using namespace constants;
        const RealVector OX {m_grid.getCellCenter(j)};

        const Real alpha = atan(-2.0);
        const Real gamma = m_eos.computeAdiabaticIndex();


        ConsState u_j;
        if (tan(alpha)*(OX[IY]-0.5)>(OX[IX]-0.5))
        {
          u_j.m(IX) = sin(alpha)*0.0 + cos(alpha)*10.0;
          u_j.m(IY) =  cos(alpha)*0.0 - sin(alpha)*10.0;
          u_j.B(IX) = sin(alpha)*5./sqrt(4.*pi) + cos(alpha)*5./sqrt(4.*pi);
          u_j.B(IY) =  cos(alpha)*5./sqrt(4.*pi) - sin(alpha)*5./sqrt(4.*pi);
          u_j.e = 20.0/(gamma-1.) + 0.5*(pow(u_j.m(IX),2)+pow(u_j.m(IY),2)+pow(u_j.B(IX),2)+pow(u_j.B(IY),2));
        }
        else
        {
          u_j.m(IX) = sin(alpha)*0.0 - cos(alpha)*10.0;
          u_j.m(IY) =  cos(alpha)*0.0 + sin(alpha)*10.0;
          u_j.B(IX) = sin(alpha)*5./sqrt(4.*pi) + cos(alpha)*5./sqrt(4.*pi);
          u_j.B(IY) =  cos(alpha)*5./sqrt(4.*pi) - sin(alpha)*5./sqrt(4.*pi);
          u_j.e = 1./(gamma-1.) + 0.5*(pow(u_j.m(IX),2)+pow(u_j.m(IY),2)+pow(u_j.B(IX),2)+pow(u_j.B(IY),2));

        }

        u_j.d = 1.0;
        u_j.m(IZ) = 0.0;
        u_j.B(IZ) = 0.0;
        
        u_j.dX=0.0;

        set(m_u, j, u_j);
    }

    Params m_params;
    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
    Real x_mid, y_mid;
};

}}
