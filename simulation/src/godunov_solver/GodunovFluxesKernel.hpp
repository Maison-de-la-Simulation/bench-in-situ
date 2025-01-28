#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "ars/HydroRiemannAllRegime.hpp"
#include "ars/MHDRiemann_AllRegime_3W.hpp"
#include "ars/MHDRiemann_AllRegime_5W.hpp"
#include "ars/MHDRiemann_ARWB_5W.hpp"
#include "ars/MHD3W_optimized.hpp"
#include "ars/HydroRiemannAllRegime_experiment_entropy.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

template < class RiemannSolver >
class FluxesAndUpdateKernel : public BaseKernel
{
public:
    using Super = BaseKernel;

    using MHD = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState = typename MHD::ConsState;
    using PrimState = typename MHD::PrimState;
    using RealVector = RealVector3d;

    using IntVector = IntVectorNd< three_d >;
    using TeamPolicy = Kokkos::TeamPolicy< Kokkos::IndexType< Int > >;
    using Team = typename TeamPolicy::member_type;
    static constexpr int ghostDepth = 0;

    FluxesAndUpdateKernel( const Params& params, const UniformGrid& grid,
                           const Array3d& u, const Array3d& q, const Kokkos::Array<Array3d, 2*three_d>& qr, Real dt )
        : m_u( u )
        , m_cu( u )
        , m_q( q )
        , m_qr( qr )
        , m_grid( grid )
        , m_eos( params.thermo )
        , m_riemann_solver( params, m_eos, m_grid )
        , m_dt( dt )
        , m_powell_st_when_low_plasma_beta (params.run.powell_st_when_low_plasma_beta)
        , m_beta_threshold(params.run.beta_threshold)
        , m_muscl_enabled(params.run.muscl_enabled)
#if defined(MPI_SESSION)
        , m_mpi_z_rank(grid.comm.getCoords(grid.comm.rank())[IZ])
#else
        , m_mpi_z_rank(0)
#endif
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const Int ix, const Int iy, const Int iz ) const
    {
      const Int j = m_grid.coordToIndex({ix,iy,iz});

      ConsState u_j = Super::getCons( m_cu, j );
      //PrimState q_j = Super::getPrim( m_q, j );

      bool st_powell = false;
      // Real beta = constants::one;
      //
      // if (m_powell_st_when_low_plasma_beta)
      // {
      //   beta = q_j.p/MHD::computeMagneticEnergy(q_j);
      //
      //   if (beta < m_beta_threshold)
      //   {
      //     st_powell = true;
      //   }
      // }

      updateAlongFace( Tags::DirX, Tags::SideL, j, u_j, st_powell );
      updateAlongFace( Tags::DirX, Tags::SideR, j, u_j, st_powell );
      updateAlongFace( Tags::DirY, Tags::SideL, j, u_j, st_powell );
      updateAlongFace( Tags::DirY, Tags::SideR, j, u_j, st_powell );
      updateAlongFace( Tags::DirZ, Tags::SideL, j, u_j, st_powell );
      updateAlongFace( Tags::DirZ, Tags::SideR, j, u_j, st_powell );

      // if ( m_powell_st_when_low_plasma_beta && (beta<m_beta_threshold) && m_muscl_enabled)
      // {
      //   Vector<three_d, Real> B_L;
      //   B_L(IX) = m_qr[0+2*IX](j, VP::IBx);
      //   B_L(IY) = m_qr[0+2*IY](j, VP::IBy);
      //   B_L(IZ) = m_qr[0+2*IZ](j, VP::IBz);
      //
      //   Vector<three_d, Real> B_R;
      //   B_R(IX) = m_qr[1+2*IX](j, VP::IBx);
      //   B_R(IY) = m_qr[1+2*IY](j, VP::IBy);
      //   B_R(IZ) = m_qr[1+2*IZ](j, VP::IBz);
      //
      //   u_j.B -= (B_R(IX) - B_L(IX))*m_dt*m_grid.m_invdl[IX]*q_j.v;
      //   u_j.B -= (B_R(IY) - B_L(IY))*m_dt*m_grid.m_invdl[IY]*q_j.v;
      //   u_j.B -= (B_R(IZ) - B_L(IZ))*m_dt*m_grid.m_invdl[IZ]*q_j.v;
      // }

      Super::set( m_u, j, u_j );
    };


    template < int idir, int iside >
    KOKKOS_FORCEINLINE_FUNCTION void
    updateAlongFace( std::integral_constant< int, idir > idir_tag,
                     std::integral_constant< int, iside > iside_tag, Int j, ConsState& u_j, bool st_powell ) const
    {
        using namespace constants;

        const Int j_L = iside == 0 ? j - m_grid.m_strides[ idir ] : j;
        const Int j_R = iside == 0 ? j : j + m_grid.m_strides[ idir ];

        /* const IntVector coord = m_grid.indexToCoord(j);
        const bool BCbot = (m_mpi_z_rank==0)&&(idir==IZ)&&(coord[IZ]==m_grid.m_ghostWidths[ IZ ])&&(iside==0);
        const bool BCtop = (m_mpi_z_rank==m_grid.m_dom[IZ]-1)&&(idir==IZ)&&(coord[IZ]==m_grid.m_nbCells[ IZ ]+m_grid.m_ghostWidths[ IZ ]-1)&&(iside==1); */

        PrimState q_L, q_R;
        //if (m_muscl_enabled && (BCtop==false) && (BCbot==false))
        if (m_muscl_enabled)
        {
            q_L = Super::getPrim( m_qr[1+2*idir], j_L );
            q_R = Super::getPrim( m_qr[0+2*idir], j_R );
        }
        else
        {
            q_L = Super::getPrim( m_q, j_L );
            q_R = Super::getPrim( m_q, j_R );
        }

        const Real side = iside == 0 ? -one : +one;

        const ConsState flux = m_riemann_solver( q_L, q_R, side, idir, st_powell );

        const Real sdtdSdV = side * m_dt * m_grid.getdSdV( j, idir_tag, iside_tag );
        u_j.d -= sdtdSdV * flux.d;
        u_j.m -= sdtdSdV * flux.m;
        u_j.e -= sdtdSdV * flux.e;
        u_j.B -= sdtdSdV * flux.B;
        u_j.dX -= sdtdSdV * flux.dX;
    }

    Array3d m_u;
    ConstArray3d m_cu;
    ConstArray3d m_q;
    Kokkos::Array<Array3d, 2*three_d> m_qr;
    UniformGrid m_grid;
    EquationOfState m_eos;
    RiemannSolver m_riemann_solver;
    Real m_dt;
    bool m_powell_st_when_low_plasma_beta;
    Real m_beta_threshold;
    bool m_muscl_enabled;
    int m_mpi_z_rank;
};

} // namespace hydro
