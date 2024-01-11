#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

class ConvectionSourceTermKernel : public BaseKernel
{
public:
    using Super = BaseKernel;

    using MHD = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState = typename MHD::ConsState;
    using PrimState = typename MHD::PrimState;

    using TeamPolicy = Kokkos::TeamPolicy< Kokkos::IndexType< Int > >;
    using Team = typename TeamPolicy::member_type;
    static constexpr int ghostDepth = 0;

    ConvectionSourceTermKernel( const Params& params, const UniformGrid& grid,
                                const Array3d& u, const Array3d& q, Real dt )
        : m_u( u )
        , m_q( q )
        , m_grid( grid )
        , m_dt( dt )
        , m_eos( params.thermo )

        , m_T_bottom(params.reader.GetReal("problem", "T_bottom", 0.0))
        , m_Bx_Bottom(params.reader.GetReal("problem",     "Bx_Bottom", 0.0))
        , m_Bz_Bottom(params.reader.GetReal("problem",     "Bz_Bottom", 0.0))
	      , m_By_Bottom(params.reader.GetReal("problem",     "By_Bottom", 0.0))

        , m_X_bottom(params.reader.GetReal("problem", "X_bottom", 0.0))

        , m_grad_T(params.reader.GetReal("problem", "grad_T", 0.0))
        , m_grad_X(params.reader.GetReal("problem", "grad_X", 0.0))

        , m_HT(params.reader.GetReal("problem", "HT", 0.0))
        , m_QA(params.reader.GetReal("problem", "QA", 0.0))
        , m_RX(params.reader.GetReal("problem", "RX", 0.0))

        , m_H_source_term(params.reader.GetBoolean("hydro", "H_source_term", false))
        , m_Q_source_term(params.reader.GetBoolean("hydro", "Q_source_term", false))
        , m_R_source_term(params.reader.GetBoolean("hydro", "R_source_term", false))

        , m_x_start(grid.m_ghostWidths[ IX ] - ghostDepth)
        , m_x_end(grid.m_nbCells[ IX ] + grid.m_ghostWidths[ IX ] + ghostDepth)
#if defined(MPI_SESSION)
        , m_mpi_z_rank(m_grid.comm.getCoords(m_grid.comm.rank())[IZ])
#else
        , m_mpi_z_rank(0)
#endif
    {
    }

    KOKKOS_INLINE_FUNCTION
    Real X0z(const Int iz) const noexcept
    {
        const Real z = (iz-m_grid.m_ghostWidths[IZ]+0.5)*m_grid.m_dl[IZ];
        return m_X_bottom + m_grad_X*z;
    }

    KOKKOS_INLINE_FUNCTION
    Real T0z(const Int iz) const noexcept
    {
        const Real z = (iz-m_grid.m_ghostWidths[IZ]+0.5)*m_grid.m_dl[IZ];
        return m_T_bottom + m_grad_T * z;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const Int ix, const Int iy, const Int iz ) const
    {
      const Int j = m_grid.coordToIndex({ix,iy,iz});

      const int iz_gloc = iz + m_mpi_z_rank*m_grid.m_nbCells[IZ];

      ConsState u_j = Super::getCons( m_u, j );
      PrimState q_j = MHD::conservativeToPrimitive(u_j, m_eos);

      const Real T_old = MHD::computeTemperature(q_j, m_eos);
      const Real Bx_old = q_j.B(IX);
      const Real By_old = q_j.B(IY);
      const Real Bz_old = q_j.B(IZ);
      const Real X_old  =  q_j.X;

      Real T_new   = T_old;
      if( m_H_source_term )
      {
        T_new = (T_old - m_dt*m_HT*T0z(iz_gloc))/(1.0-m_dt*m_HT);
      }

      Real Bx_new = Bx_old;
      Real By_new = By_old;
      Real Bz_new = Bz_old;
      if( m_Q_source_term )
      {
        Bx_new = (Bx_old - m_dt*m_QA*m_Bx_Bottom)/(1.0-m_dt*m_QA);
        By_new = (By_old - m_dt*m_QA*m_By_Bottom)/(1.0-m_dt*m_QA);
        Bz_new = (Bz_old - m_dt*m_QA*m_Bz_Bottom)/(1.0-m_dt*m_QA);
      }

      Real X_new  = X_old;
      if( m_R_source_term )
      {
        X_new = (X_old - m_dt*m_RX*X0z(iz_gloc))/(1.0-m_dt*m_RX);
      }

      q_j.B(IX) = Bx_new;
      q_j.B(IY) = By_new;
      q_j.B(IZ) = Bz_new;
      q_j.X     =  X_new;
      q_j.p = m_eos.computePressure(q_j.d, T_new, X_new);

      u_j = MHD::primitiveToConservative(q_j, m_eos);

      Super::set( m_u, j, u_j );
    } ;


    Array3d m_u;
    ConstArray3d m_q;
    UniformGrid m_grid;
    Real m_dt;
    EquationOfState m_eos;

    Real m_T_bottom;
    Real m_Bx_Bottom;
    Real m_Bz_Bottom;
    Real m_By_Bottom;
    Real m_X_bottom;

    Real m_grad_T;
    Real m_grad_X;

    Real m_HT;
    Real m_QA;
    Real m_RX;

    bool m_H_source_term;
    bool m_Q_source_term;
    bool m_R_source_term;

    Int m_x_start;
    Int m_x_end;
    int m_mpi_z_rank;
};

} // namespace hydro
