#pragma once

#include "HydroConstants.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "io/Writer.hpp"
#include "Utils.hpp"
#include "TimeStep.hpp"
#include "global_mean.hpp"
#include "Timer.hpp"

#include <memory>
#include <vector>

namespace hydro
{

class GodunovSolver: public Solver
{
    using Super           = Solver;

    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using VC              = typename MHD::VarCons;
    using VP              = typename MHD::VarPrim;
    static constexpr int nbvar = MHD::nbvar;

    static constexpr Int ghostWidth {2};

    using Array           = Array3d;
    using HostArray       = HostArray3d;

public:
    GodunovSolver(std::shared_ptr<Problem> problem);
    GodunovSolver(const GodunovSolver& x) = default;
    GodunovSolver(GodunovSolver&& x) = default;
#if defined(__INTEL_COMPILER) && __INTEL_COMPILER <= 1900
    ~GodunovSolver() override {};
#else
    ~GodunovSolver() override = default;
#endif

    Real computeTimeStep() final;
    void nextIteration(Real dt) final;
    void prepareNextOutput(Real& dt) final;
    void pdiExposeData() final;
    bool finished() const final;
    void printMonitoring(double t_tot) const final;
    bool shouldPrintInformation() const final;
    void printInformation(Real dt) const final;
    double memoryUsage() const final;
    void set_should_save() final;
    void set_time_limit_reached() final;
    void accumulate_compute_duration(const std::chrono::steady_clock::duration& duration) final;

private:
    void computeFluxesAndUpdate(Real dt);
    void compute_adjust_timestep(Real dt_type, Real dt, Real& delta_type);

private:
    std::shared_ptr<Problem> m_problem;
    std::shared_ptr<Params> m_params;
    const UniformGrid m_grid;
    std::shared_ptr<io::WriterBase> m_writer;
    bool m_should_save;
    bool m_time_limit_reached;

    const Array m_u;
    const Array m_q;
    const HostArray m_u_host;

    Kokkos::Array<Array, 2*three_d> m_qr;

    Int m_nStepmax;
    Real m_tEnd;
    TimeStep m_dt;
    global_mean m_means;
    Int m_nx;
    Int m_ny;
    Int m_nz;
    Int m_mz;
    Int m_my;

    PerformanceTimer performanceTimer;


    #if defined(MPI_SESSION)
      int m_mpi_z_rank = m_grid.comm.getCoords(m_grid.comm.rank())[IZ];
      int m_mpi_y_rank = m_grid.comm.getCoords(m_grid.comm.rank())[IY];
    #else
      int m_mpi_z_rank = 0 ;
      int m_mpi_y_rank = 0 ;
    #endif
    
};

}
