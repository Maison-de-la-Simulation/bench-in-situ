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
#include "io/WriterPDI.hpp"



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
    , performanceTimer()

{

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

void GodunovSolver::deep(double* b, double* a)
{
    Kokkos::fence();
    std::chrono::steady_clock::time_point m_start_io_deep = std::chrono::steady_clock::now();

            int* iter; PDI_access("iter", (void**)&iter, PDI_IN);
            int* freq; PDI_access("freq", (void**)&freq, PDI_IN);
            printf("########### %i ################\n", *iter);
            printf("########### %i ################\n", *freq);
//            PDI_release("freq");
//            PDI_release("iter");

    std::array<size_t, 2>* m_u_dim; PDI_access("m_u_kokkos_view_dimensions", (void**)&m_u_dim, PDI_IN);
    size_t* dim_ptr = m_u_dim->data();
    std::array<size_t, 2>* m_u_host_dim; PDI_access("m_u_host_kokkos_view_dimensions", (void**)&m_u_host_dim, PDI_IN);
    size_t* dim_host_ptr = m_u_host_dim->data();

    printf("--- dim_ptr[0] %zu ---\n", dim_ptr[0]);
    printf("--- dim_ptr[1] %zu ---\n", dim_ptr[1]);
    printf("--- dim_host_ptr[0] %zu ---\n", dim_host_ptr[0]);
    printf("--- dim_host_ptr[1] %zu ---\n", dim_host_ptr[1]);

    Kokkos::View<Real*, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u(a, dim_ptr[0], dim_ptr[1]);
    Kokkos::View<Real*, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u_host(b, dim_host_ptr[0], dim_host_ptr[1]);
    Kokkos::deep_copy(mm_u_host, mm_u);

    Kokkos::fence();
    performanceTimer.time_spent_in_deep_copy += (std::chrono::steady_clock::now() - m_start_io_deep);

    void* data_host = &mm_u_host;
//    PDI_expose("data_host", data_host, PDI_OUT);
//    free(data_host);

    PDI_multi_expose("data_GPU_event",
                     "iStep", iter, PDI_OUT,
                     NULL);
//                     "data_host", &data_host, PDI_OUT);
//                     "data_host", data_host->data(), PDI_OUT);
//                     "local_full_field", mm_u_host.data(), PDI_OUT);

    PDI_release("m_u_host_kokkos_view_dimensions");
    PDI_release("m_u_kokkos_view_dimensions");
            PDI_release("freq");
            PDI_release("iter");
}

extern "C"
{
    void copy_func() {
//        printf("*********** in copy_func ***********************\n");
        if (Session::isIOProc())
        {
            int* iter; PDI_access("iter", (void**)&iter, PDI_IN);
            int* freq; PDI_access("freq", (void**)&freq, PDI_IN);
            double* a; PDI_access("m_u", (void**)&a, PDI_IN);
            double* b; PDI_access("m_u_host", (void**)&b, PDI_IN);
            std::array<size_t, 2>* m_u_dim; PDI_access("m_u_kokkos_view_dimensions", (void**)&m_u_dim, PDI_IN);
            size_t* dim_ptr = m_u_dim->data();
            std::array<size_t, 2>* m_u_host_dim; PDI_access("m_u_host_kokkos_view_dimensions", (void**)&m_u_host_dim, PDI_IN);
            size_t* dim_host_ptr = m_u_host_dim->data();


            if (*iter % *freq == 0) {
                printf("*********** %i ***********************\n", *iter);
                printf("*********** %i ***********************\n\n", *freq);
                Kokkos::Profiling::pushRegion("I/O - Checkpoint");
                Print() << "===================== output at iteration = " << *iter << " time t = " << *time << std::endl;
                Kokkos::Profiling::pushRegion("I/O - Checkpoint - deep_copy");
                // deep(b, a);
                Kokkos::View<Real*, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u(a, dim_ptr[0], dim_ptr[1]);
                Kokkos::View<Real*, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u_host(b, dim_host_ptr[0], dim_host_ptr[1]);
                Kokkos::deep_copy(mm_u_host, mm_u);
//                Real* copied_ptr =  mm_u_host.data();
                double* copied_ptr =  mm_u_host.data();
                printf("--- after deep bis %i ---\n\n", *iter);

    // std::array<int, 3> pdi_ncells;
    // pdi_ncells[IX] = m_grid.m_nbCells[IX] * m_grid.m_dom[IX];
    // pdi_ncells[IY] = m_grid.m_nbCells[IY] * m_grid.m_dom[IY];
    // pdi_ncells[IZ] = m_grid.m_nbCells[IZ] * m_grid.m_dom[IZ];
                
    // Int& outputId = io::WriterBase::m_outputId;

    char *prefix_c_str;
    PDI_access("prefix", (void **)&prefix_c_str, PDI_IN);
    std::string prefix(prefix_c_str);
    PDI_release("prefix");

    // std::string filename = io::getFilename(prefix, outputId);
    std::string filename = io::getFilename(prefix, *iter);
    int filename_size = filename.size();

                Kokkos::Profiling::popRegion();
                Kokkos::Profiling::pushRegion("I/O - Checkpoint - write");
//                wr(b);
//                m_writer->write(b, m_grid, iter, time,
//                                m_params->thermo.gamma, m_params->thermo.mmw);
                PDI_multi_expose("data_HOST",
                                "iStep", iter, PDI_OUT,
                                "u_host", copied_ptr, PDI_OUT,
                                "m_u_host_kokkos_view_dimensions", dim_host_ptr, PDI_OUT,
                                "filename_size", &filename_size, PDI_OUT,
                                "filename", filename.data(), PDI_OUT,
                                // "grid_size", pdi_ncells.data(), PDI_OUT,
                                // "gamma", m_params->thermo.gamma, PDI_OUT,
                                NULL);
//                                "u_host", mm_u_host.data(), PDI_OUT,
//                                "u_host", (void*)(mm_u_host.data()), PDI_OUT,
//                                "m_u_host_kokkos_view_dimensions", (void*)&m_u_host_dim, PDI_OUT,
                                // "m_u_host", (void*)(m_u_host.data()), PDI_OUT,
//                                "m_u_host_kokkos_view_dimensions", (void*)&m_u_host_kokkos_view_dimensions, PDI_OUT,
                                // "data_host", &data_host, PDI_OUT);
                                // "data_host", data_host->data(), PDI_OUT);
                                // "local_full_field", mm_u_host.data(), PDI_OUT);
                Kokkos::Profiling::popRegion();
                Kokkos::Profiling::popRegion();
            }
            PDI_release("m_u_host_kokkos_view_dimensions");
            PDI_release("m_u_kokkos_view_dimensions");
            PDI_release("m_u_host");
            PDI_release("m_u");
            PDI_release("freq");
            PDI_release("iter");
        }
    }
}

void GodunovSolver::pdiExposeData()
{
  Kokkos::fence();
  std::chrono::steady_clock::time_point m_start_io = std::chrono::steady_clock::now();

#if defined(Euler_ENABLE_PDI)
    // PDI_multi_expose("data_on_GPU",
    //                  "iStep", (void*)&(Super::m_iteration), PDI_OUT,
    //                  "time", (void*)&(m_t), PDI_OUT,
    //                  NULL);

   std::array<size_t, 2> m_u_kokkos_view_dimensions = { m_u.extent(0), m_u.extent(1) };
   std::array<size_t, 2> m_u_host_kokkos_view_dimensions = { m_u_host.extent(0), m_u_host.extent(1) };

//    printf("--- pdiExposeData ---\n");
//    printf("--- mu0 %zu ---\n", m_u_kokkos_view_dimensions[0]);
//    printf("--- mu1 %zu ---\n", m_u_kokkos_view_dimensions[1]);
//    printf("--- muhost0 %zu ---\n", m_u_host_kokkos_view_dimensions[0]);
//    printf("--- muhost1 %zu ---\n\n", m_u_host_kokkos_view_dimensions[1]);
//    printf("--- &mu %p ---\n", &m_u_kokkos_view_dimensions);
//    printf("--- &muhost %p ---\n", &m_u_host_kokkos_view_dimensions);

//   PDI_multi_expose("data_on_GPU",
   PDI_multi_expose("data_GPU_event",
                    "iStep", (void*)&(Super::m_iteration), PDI_OUT,
                    "time", (void*)&(m_t), PDI_OUT,
                    "m_u", (void*)(m_u.data()), PDI_OUT,
                    "m_u_host", (void*)(m_u_host.data()), PDI_OUT,
                    "m_u_kokkos_view_dimensions", (void*)&m_u_kokkos_view_dimensions, PDI_OUT,
                    "m_u_host_kokkos_view_dimensions", (void*)&m_u_host_kokkos_view_dimensions, PDI_OUT,
                    "gamma", m_params->thermo.gamma, PDI_OUT,
                    NULL);

#endif

//     if (m_should_save)
//     {
//        Kokkos::Profiling::pushRegion("I/O - Checkpoint");
//        if(Super::m_iteration%100 == 0) Print() << "===================== output at iteration = " << Super::m_iteration << " time t = "<<Super::m_t<< std::endl;
//        Kokkos::Profiling::pushRegion("I/O - Checkpoint - deep_copy");
//        Kokkos::deep_copy(m_u_host, m_u);
//        Kokkos::Profiling::popRegion();
//        Kokkos::Profiling::pushRegion("I/O - Checkpoint - write");
//        m_writer->write(m_u_host, m_grid, Super::m_iteration, Super::m_t,
//                        m_params->thermo.gamma, m_params->thermo.mmw);
//        Kokkos::Profiling::popRegion();
//        Kokkos::Profiling::popRegion();

//    }

  Kokkos::fence();
  performanceTimer.time_spent_in_io += (std::chrono::steady_clock::now() - m_start_io);
}



bool GodunovSolver::finished() const
{
    return (Super::m_t >= m_tEnd || Super::m_iteration >= m_nStepmax || m_time_limit_reached);
}


void GodunovSolver::printMonitoring(double t_tot) const
{
    const double w_perf {Session::getNProc() * static_cast<double>(Super::m_iteration) * static_cast<double>(m_grid.nbCells()) / t_tot * 1.0E-6};
    Print() << "[RESULT] Performance " << w_perf << " Mcell-updates/s" << std::endl;
    Print() << "[RESULT] Wall_time " << t_tot << " s" << std::endl;
    Print() << performanceTimer << std::endl;
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

void GodunovSolver::accumulate_compute_duration(const std::chrono::steady_clock::duration& duration) {
  performanceTimer.time_spent_in_compute += duration;
}

}
