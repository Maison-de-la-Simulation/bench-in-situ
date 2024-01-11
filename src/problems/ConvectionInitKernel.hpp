#pragma once

#include "ConvectionParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroBaseKernel.hpp"
#include "MHDSystem.hpp"
#include <cstdlib>

// kokkos random numbers
#include <Kokkos_Random.hpp>


namespace hydro { namespace problems
{

class ConvectionInitKernel : BaseKernel
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
    using device = Kokkos::DefaultExecutionSpace;
    using GeneratorPool = Kokkos::Random_XorShift64_Pool<device>;
    using GeneratorType = GeneratorPool::generator_type;


public:
    ConvectionInitKernel(const Params& params,
                                         const ConvectionParams& prob_params,
                                         Array u,
                                         const UniformGrid& grid)
        : m_params {params}
        , m_prob_params {prob_params}
        , m_eos {params.thermo}
        , m_u {u}
        , m_grid {grid}
        , m_kB {params.thermo.kB}
        , m_random_perturbation {params.run.random_perturbation}
        , m_three_d_perturbation {params.run.three_d_perturbation}
        #if defined(MPI_SESSION)
        , m_mpi_z_rank {m_grid.comm.getCoords(m_grid.comm.rank())[IZ]}
        , m_mpi_rank {m_grid.comm.rank()}
        #else
        , m_mpi_z_rank {0}
        , m_mpi_rank {0}
        #endif
        , m_random_offset{static_cast<uint64_t>(m_mpi_rank+50)}
        , rand_pool{Kokkos::Random_XorShift64_Pool<device>( m_random_offset ) }

    {
    }

    static void apply(const Params& params,
                      const ConvectionParams& prob_params,
                      Array u,
                      const UniformGrid& grid)
    {
      ConvectionInitKernel kernel {params, prob_params, u, grid};

      Int start_x = 0;
      Int start_y = 0;

      Int end_x = grid.m_nbCells[IX]+2*grid.m_ghostWidths[IX];
      Int end_y = grid.m_nbCells[IY]+2*grid.m_ghostWidths[IY];

      Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy ({start_x, start_y}, {end_x, end_y});
      Kokkos::parallel_for(policy, kernel);
    }

    KOKKOS_INLINE_FUNCTION
    Real u0(Int index) const
    {
        const RealVector OX{m_grid.getCellCenter(index)};
        Real u0;
        u0 = m_prob_params.u0_bottom + m_prob_params.grad_u0 * OX[IZ];
        return u0;
    }

    KOKKOS_INLINE_FUNCTION
    Real T0(Int index) const
    {
        const RealVector OX{m_grid.getCellCenter(index)};
        Real T;

          T = m_prob_params.T_bottom + m_prob_params.grad_T * (OX[IZ]);

        return T;
    }

    KOKKOS_INLINE_FUNCTION
    Real X0(Int index) const
    {
      const RealVector OX{m_grid.getCellCenter(index)};
      Real X;

        X = m_prob_params.X_bottom + m_prob_params.grad_X * (OX[IZ]);

      return X;
    }

    KOKKOS_INLINE_FUNCTION
    Real X0z(Int iz) const
    {
      Real X;

        Real z = (iz-m_grid.m_ghostWidths[IZ]+0.5)*m_grid.m_dl[IZ];
        X = m_prob_params.X_bottom + m_prob_params.grad_X * z;

      return X;
    }

    KOKKOS_INLINE_FUNCTION
    Real T0z(Int iz) const
    {
        Real T;

          Real z = (iz-m_grid.m_ghostWidths[IZ]+0.5)*m_grid.m_dl[IZ];
          T = m_prob_params.T_bottom + m_prob_params.grad_T * z;

        return T;
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
    void operator()(Int ix, Int iy) const
    {
        using namespace constants;

            for (Int iz=m_grid.m_ghostWidths[IZ];
                 iz<m_grid.m_nbCells[IZ]+m_grid.m_ghostWidths[IZ]; ++iz)
            {

                const IntVector coord {{ix, iy, iz}};
                const Int j {m_grid.coordToIndex(coord)};
                const IntVector coordm1 {{ix, iy, iz-1}};
                const Int jm1 {m_grid.coordToIndex(coordm1)};

                PrimState q_j;
                Real rho_0 {m_prob_params.density_bottom};

                int iz_gloc = iz + m_mpi_z_rank*m_grid.m_nbCells[IZ];

                const Real alpha = 0.5*(phi(j) - phi(jm1))/m_kB;

                Real rho_iz=rho_0;

                if (iz_gloc>=m_grid.m_ghostWidths[IZ])
                {
                  for(int k = m_grid.m_ghostWidths[IZ]; k<=iz_gloc; k++)
                  {
                     Real mu_km1 = m_eos.computeMeanMolecularWeight(X0z( k-1 ));
                     Real mu_k   = m_eos.computeMeanMolecularWeight(X0z( k   ));
                     Real T0_km1 =                                  T0z( k-1 );
                     Real T0_k =                                    T0z( k   );

                     rho_iz*=( T0_km1/mu_km1 - alpha)/(T0_k/mu_k + alpha);
                  }
                }

                q_j.d = rho_iz;

                q_j.p = m_eos.computePressure(q_j.d, T0z(iz_gloc), X0z(iz_gloc));
                q_j.v( IX ) = zero;
                q_j.v( IY ) = zero;
                q_j.v( IZ ) = zero;
                q_j.B( IX ) = m_prob_params.Bx_bottom;
                q_j.B( IY ) = m_prob_params.By_bottom;
                q_j.B( IZ ) = m_prob_params.Bz_bottom;
                q_j.X = X0z(iz_gloc);


                const RealVector OX{m_grid.getCellCenter(j)};
                const Real Kx = constants::pi*m_prob_params.kx;
                const Real Ky = constants::pi*m_prob_params.ky;
                const Real Kz = constants::pi*m_prob_params.kz;

                // get random number state
                if (m_random_perturbation)
                {
                  GeneratorType rand_gen = rand_pool.get_state();
                  q_j.v( IZ ) =  m_prob_params.amplitude_seed*rand_gen.drand()*exp(-(OX[IZ]-0.5)*(OX[IZ]-0.5)*40);
                  rand_pool.free_state(rand_gen);
                }
                else
                {
                  if (m_three_d_perturbation)
                  {
                    q_j.v( IZ ) =  m_prob_params.amplitude_seed*sin(Kx*OX[IX])*sin(Kz*OX[IZ])*sin(Ky*OX[IY]);
                  }
                  else
                  {
                    q_j.v( IZ ) =  m_prob_params.amplitude_seed*sin(Kx*OX[IX])*sin(Kz*OX[IZ]);
                  }
                }



                const ConsState u_j {MHD::primitiveToConservative(q_j, m_eos)};

                set(m_u, j, u_j);

              // if ((iz == m_grid.m_nbCells[IZ]+m_grid.m_ghostWidths[IZ]-1)&&(ix==1)&&(iy==0))
              // {
              //
              //   const IntVector coord_top {{ix, iy, iz}};
              //   const Int j_top {m_grid.coordToIndex(coord_top)};
              //
              //   const IntVector coord_bot {{ix, iy, m_grid.m_ghostWidths[IZ]}};
              //   const Int j_bot {m_grid.coordToIndex(coord_bot)};
              //
              //   PrimState q_top = MHD::conservativeToPrimitive(Super::getCons(m_u, j_top), m_eos);
              //   PrimState q_bot = MHD::conservativeToPrimitive(Super::getCons(m_u, j_bot), m_eos);
              //
              //   Real zmax= 1.0;
              //
              //   Real log_p_bot = log(q_bot.p);
              //   Real log_p_top = log(q_top.p);
              //
              //   Real log_T_bot = log(MHD::computeTemperature(q_bot, m_eos));
              //   Real log_T_top = log(MHD::computeTemperature(q_top, m_eos));
              //
              //   Real log_mu_bot = log(m_eos.computeMeanMolecularWeight(q_bot.X));
              //   Real log_mu_top = log(m_eos.computeMeanMolecularWeight(q_top.X));
              //
              //   Real hp = -1.0/((log_p_top-log_p_bot)/zmax); //Only valid for hmax=0.5
              //   Real nabla_T = -hp*((log_T_top - log_T_bot)/zmax);
              //   Real nabla_ad = (m_eos.computeAdiabaticIndex() - 1.0)/m_eos.computeAdiabaticIndex();
              //
              //   Real nabla_mu = -hp*((log_mu_top - log_mu_bot)/zmax);
              //
              //   Real HT = m_prob_params.HT;
              //   Real QA = m_prob_params.QA;
              //   Real RX = m_prob_params.RX;
              //
              //   Real g_z = m_prob_params.g_z;
              //   Real rho = 0.5*(q_top.d + q_bot.d);
              //   Real B0 = m_prob_params.Bx_bottom;
              //
              //   Real kx2 = Kx*Kx;
              //   Real ky2 = Ky*Ky;
              //   Real kz2 = Kz*Kz;
              //
              //   Real k2 = kx2+ky2+kz2;
              //
              //   Real B02 = B0*B0;
              //
              //   Real CB = -(k2*kx2/(kx2+ky2))*(hp/(g_z*rho));
              //
              //   Real adiabatic = nabla_T - nabla_ad - nabla_mu - CB*B02;
              //   Real diabatic = (nabla_T - nabla_ad)*(RX+QA) - nabla_mu*(QA+HT) - CB*B02*(HT+RX);
              //   Real ddiabatic = (nabla_T - nabla_ad)*RX*QA - nabla_mu*QA*HT - CB*B02*HT*RX;
              //
              //
              //   printf(" \n HT = %.4f \n", HT);
              //   printf(" QA = %.4f \n", QA);
              //   printf(" RX = %.4f \n", RX);
              //
              //   printf(" nabla_T -nabla_ad= %.4f \n", nabla_T-nabla_ad);
              //   printf(" nabla_ad = %.4f \n", nabla_ad);
              //   printf(" nabla_mu = %.4f \n", nabla_mu);
              //   printf(" CkB02 = %.4f \n", CB*B02);
              //
              //   printf( "B target = %.6f \n", sqrt(  -nabla_mu/(CB*m_prob_params.phi_target)) );
              //   printf( "Phi = %.4f \n", -nabla_mu/(CB*B02) );
              //   printf( "r = %.4f \n", RX/QA );
              //
              //   printf(" adiabatic criterion= %.4f \n", adiabatic);
              //   printf(" diabatic criterion= %.4f \n", -diabatic);
              //   printf(" ddiabatic criterion= %.5f \n", ddiabatic);
              //
              // }
            }
    }

    Params m_params;
    ConvectionParams m_prob_params;
    EquationOfState m_eos;
    Array m_u;
    UniformGrid m_grid;
    Real m_kB;
    bool m_random_perturbation;
    bool m_three_d_perturbation;
    int m_mpi_z_rank;
    int m_mpi_rank;
    uint64_t m_random_offset;
    GeneratorPool rand_pool;


};

}}
