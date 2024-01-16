#include "Godunov.hpp"

#include "ConservativeToPrimitiveExecution.hpp"
#include "DistributedMemorySession.hpp"
#include "FluxesAndUpdateKernelDispatch.hpp"
#include "HydroConstants.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "io/Reader.hpp"
#include "io/Writer.hpp"
#include "Print.hpp"
#include "ConvectionSourceTermExecution.hpp"
#include "MagneticResistivityExecution.hpp"
#include "TimeStep.hpp"
#include "TimeStepExecution.hpp"
#include "Utils.hpp"
#include "MusclReconstructionExecution.hpp"
#include "global_meanExecution.hpp"
#include "vp2Execution.hpp"



#if defined(Euler_ENABLE_PDI)
#include<pdi.h>
#endif

#include <chrono>
#include <iomanip>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace hydro
{

GodunovSolver::GodunovSolver(std::shared_ptr<Problem> problem)
    : Solver           {}
    , m_problem        {problem}
    , m_params         {problem->m_params}
    , m_grid           {m_params->mesh.low, m_params->mesh.up, m_params->mesh.nbCells,
                        m_params->mesh.dom, ghostWidth}
    , m_writer         {}
    , m_should_save    {false}
    , m_time_limit_reached  {false}
    , m_u              ("U", m_grid.nbCells())
    , m_q              ("Q", m_grid.nbCells())
    , m_u_host         {Kokkos::create_mirror(m_u)}
    , m_qr             {}
    , m_nStepmax {m_params->run.nStepmax}
    , m_tEnd {m_params->run.tEnd}
    , m_dt {}

    , m_nx {m_params->mesh.nbCells[IX]}
    , m_ny {m_params->mesh.nbCells[IY]}
    , m_nz {m_params->mesh.nbCells[IZ]}
    , m_mz {m_params->mesh.dom[IZ]}
    , m_my {m_params->mesh.dom[IY]}

{
    //slice stuff
    m_iz_middle_gloc =m_mz*m_nz/2-1;
    m_contains_middle_z = ((m_mpi_z_rank*m_nz<=m_iz_middle_gloc)&&(m_iz_middle_gloc<=(m_mpi_z_rank+1)*m_nz-1));
    iz_middle= m_iz_middle_gloc%m_nz;

    m_iy_middle_gloc =m_my*m_ny/2-1;
    m_contains_middle_y = ((m_mpi_y_rank*m_ny<=m_iy_middle_gloc)&&(m_iy_middle_gloc<=(m_mpi_y_rank+1)*m_ny-1));
    iy_middle= m_iy_middle_gloc%m_ny;
    //--

    std::vector<std::string> var_names = MHD::cons_names();
    std::vector<std::pair<int, std::string>> variables_to_save;
    for (int ivar = 0; ivar < nbvar; ++ivar)
    {
        variables_to_save.push_back(std::make_pair(ivar, var_names[ivar]));
    }

    m_writer = io::WriterFactory::New(m_grid, *m_params,
                                           m_params->output.type,
                                           m_params->output.prefix,
                                           variables_to_save);
    if(m_params->run.muscl_enabled)
    {
      for (int idim=0; idim<three_d; ++idim)
      {
          m_qr[0+2*idim] = Array("qr", static_cast<typename Array::size_type>(m_grid.nbCells()));
          m_qr[1+2*idim] = Array("qr", static_cast<typename Array::size_type>(m_grid.nbCells()));
        }
    }

    if (m_params->run.restart)
    {
        Int outputId = -1;
        Int restartId = -1;
        io::Reader reader(m_grid, *m_params, variables_to_save);
        reader.read(m_u_host, m_grid, Super::m_iteration, Super::m_t, outputId, restartId);
        Kokkos::deep_copy(m_u, m_u_host);
        m_writer->setOutputId(++outputId);
        m_writer->setRestartId(++restartId);
        m_should_save = false;
    }
    else
    {
        m_problem->initialize(m_u, m_grid);
        m_should_save = (m_params->output.nOutput > 0 ) || (m_params->output.dt_io > 0);
    }
    m_problem->make_boundaries(m_u, m_grid);
    ExecuteConservativeToPrimitive(*m_params, m_grid, m_u, m_q);
}


Real GodunovSolver::computeTimeStep()
{
    m_dt = ExecuteTimeStep(*m_params, m_grid, m_q);
    const Real dt_s {m_dt.min()};
    const Real dt {Super::m_t + dt_s > m_tEnd ? utils::adjust(Super::m_t, m_tEnd) : dt_s};
    return dt;
}


void GodunovSolver::nextIteration(Real dt)
{
    if (m_params->hydro.hydro_enabled)
    {
        if(m_params->run.muscl_enabled)
        {
          ExecuteMusclReconstruction(*m_params, m_grid, m_q, m_qr, dt);
        }
        FluxesAndUpdateKernelDispatch(*m_params, m_grid, m_u, m_q, m_qr, dt);

    }

    if (m_params->hydro.convection_source_term_enabled)
    {
        ExecuteConvectionSourceTerm(*m_params, m_grid, m_u, m_q, dt);
    }

    if (m_params->hydro.magnetic_resistivity_enabled)
    {   
        //To update from value at time n, comment these brackets. Else, it's taken after the FV update
        {
        m_problem->make_boundaries(m_u, m_grid);
        ExecuteConservativeToPrimitive(*m_params, m_grid, m_u, m_q);
        }
        ExecuteMagneticResistivity(*m_params, m_grid, m_u, m_q, dt);
    }

    // fill ghost cell in data_in
    Kokkos::Profiling::pushRegion("Boundaries");
    m_problem->make_boundaries(m_u, m_grid);
    Kokkos::Profiling::popRegion();

    // convert conservative variable into primitives ones for the entire domain
    ExecuteConservativeToPrimitive(*m_params, m_grid, m_u, m_q);

    Super::m_t += dt;
    Super::m_iteration++;
}

void GodunovSolver::prepareNextOutput(Real& dt)
{
    m_should_save = false;
    
    auto dt_io = m_params->output.dt_io;
    
    auto delta_io=std::numeric_limits<Real>::infinity();
    
    if (dt_io > constants::zero) {compute_adjust_timestep(dt_io, dt, delta_io);}
    
    dt=std::min({dt, delta_io});

    if (dt == delta_io) {m_should_save=true;}
    
    if (m_t + dt >= m_tEnd)
    {
        m_should_save = true;
    }

    if ((m_params->output.nOutput > 0) && ((Super::m_iteration + 1) % m_params->output.nOutput == 0))
    {
        m_should_save = true;
    }
    
}

void GodunovSolver::pdiExposeData()
{
#if defined(Euler_ENABLE_PDI)
    PDI_multi_expose("data_on_GPU",
                     "iStep", (void*)&(Super::m_iteration), PDI_OUT,
                     "time", (void*)&(m_t), PDI_OUT,
                     NULL);
#endif

    if (m_should_save)
    {
        Kokkos::Profiling::pushRegion("I/O - Checkpoint");
        Print() << "===================== output at iteration = " << Super::m_iteration << " time t = "<<Super::m_t<< std::endl;
        Kokkos::Profiling::pushRegion("I/O - Checkpoint - deep_copy");
        Kokkos::deep_copy(m_u_host, m_u);
        Kokkos::Profiling::popRegion();
        Kokkos::Profiling::pushRegion("I/O - Checkpoint - write");
        m_writer->write(m_u_host, m_grid, Super::m_iteration, Super::m_t,
                        m_params->thermo.gamma, m_params->thermo.mmw);
        Kokkos::Profiling::popRegion();
        Kokkos::Profiling::popRegion();

    }

}


bool GodunovSolver::finished() const
{
    return (Super::m_t >= m_tEnd || Super::m_iteration >= m_nStepmax || m_time_limit_reached);
}


void GodunovSolver::printMonitoring(double t_tot) const
{
    const double w_perf {Session::getNProc() * static_cast<double>(Super::m_iteration) * static_cast<double>(m_grid.nbCells()) / t_tot * 1.0E-6};
    Print() << "Perf (Wall clock): " << w_perf << " Mcell-updates/s\n";
}


bool GodunovSolver::shouldPrintInformation() const
{
    if (m_params->run.info == 0)
    {
        return false;
    }

    if (finished())
    {
        return true;
    }

    if (m_should_save)
    {
        return true;
    }

    return Super::m_iteration%m_params->run.info == 0;
}


void GodunovSolver::printInformation(Real dt) const
{
    Print oss {};
    oss << std::setprecision(std::numeric_limits<Real>::digits10);
    oss << std::scientific;
    oss << "Step n=";
    oss << std::setw(std::numeric_limits<int>::digits10) << std::setfill('.') << Super::m_iteration;
    oss << "; time t=" << Super::m_t;
    oss << " [" << std::setw(5) << std::setfill(' ') << std::setprecision(1) << std::fixed << 100.0*static_cast<double>(Super::m_t/m_tEnd) << "%]\n";
    oss << std::setprecision(std::numeric_limits<Real>::digits10);
    oss << std::scientific;
    oss << " * Time step CFL dt=" << m_dt.min() <<" Time step used dt="<<dt<< std::endl;
}


double GodunovSolver::memoryUsage() const
{
    auto memory = m_u.span() + m_q.span();
    for (int idim=0; idim<three_d; ++idim)

    {
      memory += m_qr[0+2*idim].span();
      memory += m_qr[1+2*idim].span();
    }

    return static_cast<double>(memory * sizeof(Real));
}

void GodunovSolver::set_should_save()
{
  m_should_save=true;
}

void GodunovSolver::set_time_limit_reached()
{
  m_time_limit_reached=true;
}

void GodunovSolver::compute_adjust_timestep(Real dt_type, Real dt, Real& delta_type)
{
 
    // Next physical time to do output
    auto t_type = (std::floor(Super::m_t / dt_type) + constants::one)*dt_type;
    if (Super::m_t + dt >= t_type)
    {
        delta_type = utils::adjust(Super::m_t, t_type,
        [dt_type](Real v1, Real v2, Real delta)
        {
            return (delta + v1 < v2)
                || ((std::floor((v1+delta) / dt_type) + constants::one)*dt_type) <= v2;
        });


        if (dt < delta_type)
        {
         throw std::runtime_error("Time step is increasing whereas it should decrease.\n");
        }
    } 
}

}
