#pragma once

#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroBaseKernel.hpp"
#include "MHDSystem.hpp"
#include "ConvectionParams.hpp"

namespace hydro { namespace problems
{

class ConvectionBoundariesKernel : public BoundariesKernel
{
    using Super           = BaseKernel;

    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState       = typename MHD::ConsState;
    using PrimState       = typename MHD::PrimState;
    using RealVector      = RealVector3d;
    using VC              = typename MHD::VarCons;
    using VP              = typename MHD::VarPrim;
    static constexpr int nbvar = MHD::nbvar;

    using Array           = Array3d;
    using IntVector       = IntVectorNd<three_d>;

public:
    struct downTag {};
    struct upTag {};

    ConvectionBoundariesKernel(Array u, const UniformGrid& grid,
                                               const Params& params,
                                               const ConvectionParams& prob_params,
                                               Int side)
        : BoundariesKernel {u, grid, 0, side}
        , m_u {u}
        , m_grid {grid}
        , m_params {params}
        , m_prob_params {prob_params}
        , m_eos {params.thermo}
        , m_kB {params.thermo.kB}{}


    static void apply(Array u, const UniformGrid& grid,
                      const Params& params,
                      const ConvectionParams& prob_params,
                      Int side)
    {
        ConvectionBoundariesKernel kernel {u, grid, params, prob_params, side};

        Int start_x = 0;
        Int start_y = 0;

        Int end_x = grid.m_nbCells[IX]+2*grid.m_ghostWidths[IX];
        Int end_y = grid.m_nbCells[IY]+2*grid.m_ghostWidths[IY];

        if (side == 0)
        {
          Kokkos::MDRangePolicy<Kokkos::Rank<2>, downTag> policy ({start_x, start_y}, {end_x, end_y});
          Kokkos::parallel_for(policy, kernel);
        }
        else
        {
          Kokkos::MDRangePolicy<Kokkos::Rank<2>, upTag > policy ({start_x, start_y}, {end_x, end_y});
          Kokkos::parallel_for(policy, kernel);
        }
    }

    KOKKOS_INLINE_FUNCTION
    Real phi(Int index) const
    {
        const RealVector OX{m_grid.getCellCenter(index)};
        Real phi;
        phi = -m_params.hydro.g[IZ] * OX[IZ];
        return phi;
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const downTag&, Int ix , Int iy) const
    {
        using namespace constants;

          for (Int iz {m_grid.m_ghostWidths[IZ]-1}; iz>=0; --iz)
          {
              const IntVector coord {{ix, iy, iz}};
              const Int j {m_grid.coordToIndex(coord)};
              const IntVector coord_k {{ix, iy, iz+1}};
              const Int k {m_grid.coordToIndex(coord_k)};
              const IntVector coord_kk {{ix, iy, iz+2}};
              const Int kk {m_grid.coordToIndex(coord_kk)};

              const ConsState u_k {getCons(m_u, k)};
              const PrimState q_k {MHD::conservativeToPrimitive(u_k, m_eos)};
              const ConsState u_kk {getCons(m_u, kk)};
              const PrimState q_kk {MHD::conservativeToPrimitive(u_kk, m_eos)};

              const Real T_k   = MHD::computeTemperature( q_k, m_eos);
              const Real T_kk  = MHD::computeTemperature(q_kk, m_eos);

              const Real T_j = 2*T_k     - T_kk  ;
              //const Real X_j = 2*q_k.X - 2*q_kk.X;

              Real X_j = 2*q_k.X - q_kk.X;

              X_j = X_j < 0.0 ? 0.0 : X_j ;
              X_j = X_j > 1.0 ? 1.0 : X_j ;

              //const Real X_j = q_k.X;

              const Real mu_k = m_eos.computeMeanMolecularWeight(q_k.X);
              const Real mu_j = m_eos.computeMeanMolecularWeight(X_j);

              PrimState q_j;

              const Real alpha = 0.5*(phi(k) - phi(j))/m_kB;

              q_j.d = q_k.d * (T_k/mu_k + alpha)/(T_j/mu_j - alpha);
              q_j.p = m_eos.computePressure(q_j.d, T_j, X_j);

              q_j.v( IX ) = -q_k.v( IX );
              q_j.v( IY ) = -q_k.v( IY );
              q_j.v( IZ ) = -q_k.v( IZ );

              q_j.B( IX ) =  q_k.B( IX );
              q_j.B( IY ) =  q_k.B( IY );
              q_j.B( IZ ) =  q_k.B( IZ );

              q_j.X= X_j;

              const ConsState u_j {MHD::primitiveToConservative(q_j, m_eos)};

              set(m_u, j, u_j);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const upTag&, Int ix, Int iy) const
    {
        using namespace constants;

        for (Int iz=m_grid.m_nbCells[IZ]+m_grid.m_ghostWidths[IZ];
             iz<m_grid.m_nbCells[IZ]+2*m_grid.m_ghostWidths[IZ]; ++iz)
        {
            const IntVector coord {{ix, iy, iz}};
            const Int j {m_grid.coordToIndex(coord)};
            const IntVector coord_k {{ix, iy, iz-1}};
            const Int k {m_grid.coordToIndex(coord_k)};
            const IntVector coord_kk {{ix, iy, iz-2}};
            const Int kk {m_grid.coordToIndex(coord_kk)};

            const ConsState u_k {getCons(m_u, k)};
            const PrimState q_k {MHD::conservativeToPrimitive(u_k, m_eos)};
            const ConsState u_kk {getCons(m_u, kk)};
            const PrimState q_kk {MHD::conservativeToPrimitive(u_kk, m_eos)};

            const Real T_k   = MHD::computeTemperature( q_k, m_eos);
            const Real T_kk  = MHD::computeTemperature(q_kk, m_eos);

            const Real T_j = 2*T_k - T_kk  ;

            Real X_j = 2*q_k.X - q_kk.X;

            X_j = X_j < 0.0 ? 0.0 : X_j ;
            X_j = X_j > 1.0 ? 1.0 : X_j ;

            //const Real X_j = q_k.X;

            const Real mu_k = m_eos.computeMeanMolecularWeight(q_k.X);
            const Real mu_j = m_eos.computeMeanMolecularWeight(X_j);

            PrimState q_j;

            const Real alpha = 0.5*(phi(j) - phi(k))/m_kB;

            q_j.d = q_k.d * (T_k/mu_k - alpha)/(T_j/mu_j + alpha);
            q_j.p = m_eos.computePressure(q_j.d, T_j, X_j);

            q_j.v( IX ) = -q_k.v( IX );
            q_j.v( IY ) = -q_k.v( IY );
            q_j.v( IZ ) = -q_k.v( IZ );

            q_j.B( IX ) =  q_k.B( IX );
            q_j.B( IY ) =  q_k.B( IY );
            q_j.B( IZ ) =  q_k.B( IZ );

            q_j.X = X_j;

            const ConsState u_j {MHD::primitiveToConservative(q_j, m_eos)};

            set(m_u, j, u_j);

        }
    }

    Array m_u;
    UniformGrid m_grid;
    Params m_params;
    ConvectionParams m_prob_params;
    EquationOfState m_eos;
    Real m_kB;
};

}}
