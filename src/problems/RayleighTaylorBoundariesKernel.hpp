#pragma once

#include "HydroUniformGrid.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "RayleighTaylorParams.hpp"

namespace hydro { namespace problems
{

class RayleighTaylorBoundariesKernel : public BoundariesKernel
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
    struct downTag {};
    struct upTag {};

    RayleighTaylorBoundariesKernel(Array u, const UniformGrid& grid, const Params& params,
                                   Int idim, Int side)
        : BoundariesKernel {u, grid, idim, side}
        , m_grid {grid}
        , m_params {params}
        , m_eos {params.thermo}
    {
    }

    static void apply(Array u, const UniformGrid& grid, const Params& params,
                      Int idim, Int side)
    {
        Int start_x = 0;
        Int start_z = 0;

        Int end_x = grid.m_nbCells[IX]+2*grid.m_ghostWidths[IX];
        Int end_z = grid.m_nbCells[IZ]+2*grid.m_ghostWidths[IZ];

        if (side==0)
        {
            Kokkos::MDRangePolicy<Kokkos::Rank<2>, downTag> policy ({start_x, start_z},{end_x, end_z});
            RayleighTaylorBoundariesKernel kernel {u, grid, params, idim, side};
            Kokkos::parallel_for(policy, kernel);
        }
        else
        {
            Kokkos::MDRangePolicy<Kokkos::Rank<2>, upTag> policy ({start_x, start_z},{end_x, end_z});
            RayleighTaylorBoundariesKernel kernel {u, grid, params, idim, side};
            Kokkos::parallel_for(policy, kernel);
        }
    }

    KOKKOS_INLINE_FUNCTION
    Real phi(Int index) const
    {
        const RealVector OX{m_grid.getCellCenter(index)};
        Real phi{};
        for (Int idim2 = 0; idim2 < two_d; ++idim2)
        {
            phi -= m_params.hydro.g[idim2] * OX[idim2];
        }
        return phi;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(downTag, Int ix, Int iz) const
    {
        using namespace constants;

        for (Int iy_=m_grid.m_ghostWidths[IY]; iy_>=1; --iy_)
        {
            const Int iy {iy_-1};
            const IntVector coord {{ix, iy, iz}};
            const Int j {m_grid.coordToIndex(coord)};
            const IntVector coord0 {{ix, iy+1, iz}};
            const Int k {m_grid.coordToIndex(coord0)};
            const ConsState u_k {getCons(m_u, k)};
            const PrimState q_k {Euler::conservativeToPrimitive(u_k, m_eos)};

            PrimState q_j;
            q_j.d = q_k.d;
            q_j.p = q_k.p + half*(q_k.d + q_j.d)*(phi(k) - phi(j));
            q_j.v( IX ) = q_k.v( IX );
            q_j.v( IY ) = - q_k.v( IY );
            q_j.v( IZ ) = - 0.0;
            q_j.B( IX ) = 0.0;
            q_j.B( IY ) = 0.0;
            q_j.B( IZ ) = 0.0;
            q_j.X=0.0;


            const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
            set(m_u, j, u_j);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(upTag, Int ix, Int iz) const
    {
        using namespace constants;

        for (Int iy=m_grid.m_nbCells[IY]+m_grid.m_ghostWidths[IY];
             iy<m_grid.m_nbCells[IY]+2*m_grid.m_ghostWidths[IY]; ++iy)
        {
            const IntVector coord {{ix, iy, iz}};
            const Int j {m_grid.coordToIndex(coord)};
            const IntVector coord0 {{ix, iy-1, iz}};
            const Int k {m_grid.coordToIndex(coord0)};
            const ConsState u_k {getCons(m_u, k)};
            const PrimState q_k {Euler::conservativeToPrimitive(u_k, m_eos)};

            PrimState q_j;
            q_j.d = q_k.d;
            q_j.p = q_k.p + half*(q_k.d + q_j.d)*(phi(k) - phi(j));
            q_j.v( IX ) = q_k.v( IX );
            q_j.v( IY ) = - q_k.v( IY );
            q_j.v( IZ ) = 0.0;
            q_j.B( IX ) = 0.0;
            q_j.B( IY ) = 0.0;
            q_j.B( IZ ) = 0.0;
            q_j.X=0.0;


            const ConsState u_j {Euler::primitiveToConservative(q_j, m_eos)};
            set(m_u, j, u_j);
        }
    }

    UniformGrid m_grid;
    Params m_params;
    EquationOfState m_eos;
};

}}
