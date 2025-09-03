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
#include "io/WriterGpuPDI.hpp"



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


extern "C"
{
    void copy_func() {
        int* iter; PDI_access("iter", (void**)&iter, PDI_IN);
        int* freq; PDI_access("freq", (void**)&freq, PDI_IN);
        Real* time; PDI_access("time", (void**)&time, PDI_IN);
        Real* a; PDI_access("m_u", (void**)&a, PDI_IN); //Real, and not just double
        Real* b; PDI_access("m_u_host", (void**)&b, PDI_IN);
        std::array<size_t, 2>* m_u_dim; PDI_access("m_u_kokkos_view_dimensions", (void**)&m_u_dim, PDI_IN);
        size_t* dim_ptr = m_u_dim->data();
        std::array<size_t, 2>* m_u_host_dim; PDI_access("m_u_host_kokkos_view_dimensions", (void**)&m_u_host_dim, PDI_IN);
        size_t* dim_host_ptr = m_u_host_dim->data();

        if (*iter % *freq == 0) {
            // printf("*********** %i ***********************\n", *iter);
            // printf("*********** %i ***********************\n\n", *freq);
            Kokkos::Profiling::pushRegion("I/O - Checkpoint");
            Print() << "===================== output at iteration = " << *iter << " time t = " << *time << std::endl;
            Kokkos::Profiling::pushRegion("I/O - Checkpoint - deep_copy");
            Kokkos::View<Real**, Kokkos::LayoutLeft> mm_u(a, dim_ptr[0], dim_ptr[1]);
            Kokkos::View<Real**, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u_host(b, dim_host_ptr[0], dim_host_ptr[1]);
            Kokkos::deep_copy(mm_u_host, mm_u);

            Real* copied_ptr = const_cast<Real*>(mm_u_host.data());
            // printf("--- after deep bis %i ---\n\n", *iter);

            char *prefix_c_str;
            PDI_access("prefix", (void **)&prefix_c_str, PDI_IN);
            // printf("prefix_c_str %s \n", prefix_c_str);
            std::string prefix(prefix_c_str);
            // printf("prefix %s \n", prefix.c_str());
            PDI_release("prefix");

            std::string filename = io::WriterGpuPDI::getFilename(prefix, *iter);
            int filename_size = filename.size();

            Kokkos::Profiling::popRegion();
            Kokkos::Profiling::pushRegion("I/O - Checkpoint - write");
//                wr(b);
//                m_writer->write(b, m_grid, iter, time,
//                                m_params->thermo.gamma, m_params->thermo.mmw);
            PDI_multi_expose("data_HOST",
                            "iStep", iter, PDI_OUT,
                            "local_full_field", copied_ptr, PDI_OUT, // u_host
                            "m_u_host_kokkos_view_dimensions", dim_host_ptr, PDI_OUT,
                            "filename_size", &filename_size, PDI_OUT,
                            "filename", filename.data(), PDI_OUT,
                            NULL);
            Kokkos::Profiling::popRegion();
            Kokkos::Profiling::popRegion();
        }
        PDI_release("m_u_host_kokkos_view_dimensions");
        PDI_release("m_u_kokkos_view_dimensions");
        PDI_release("m_u_host");
        PDI_release("m_u");
        PDI_release("time");
        PDI_release("freq");
        PDI_release("iter");
    }

    void before_func() {
        int* iter; PDI_access("iter", (void**)&iter, PDI_IN);
        int* freq; PDI_access("freq", (void**)&freq, PDI_IN);
        int* time; PDI_access("time", (void**)&time, PDI_IN);
        Real* a; PDI_access("m_u", (void**)&a, PDI_IN); //Real, and not just double
        Real* b; PDI_access("m_u_host", (void**)&b, PDI_IN);
        std::array<size_t, 2>* m_u_dim; PDI_access("m_u_kokkos_view_dimensions", (void**)&m_u_dim, PDI_IN);
        size_t* dim_ptr = m_u_dim->data();
        std::array<size_t, 2>* m_u_host_dim; PDI_access("m_u_host_kokkos_view_dimensions", (void**)&m_u_host_dim, PDI_IN);
        size_t* dim_host_ptr = m_u_host_dim->data();

        if (*iter % *freq == 0) {
            // printf("*********** %i ***********************\n", *iter);
            // printf("*********** %i ***********************\n\n", *freq);
            Kokkos::Profiling::pushRegion("I/O - Checkpoint");
            // Print() << "===================== output at iteration = " << *iter << " time t = " << *time << std::endl;
            Kokkos::Profiling::pushRegion("I/O - Checkpoint - deep_copy");
            // Kokkos::deep_copy(b, a);
            // Kokkos::View<Real**, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u(a, dim_ptr[0], dim_ptr[1]);
            Kokkos::View<Real**, Kokkos::LayoutLeft> mm_u(a, dim_ptr[0], dim_ptr[1]);
            Kokkos::View<Real**, Kokkos::LayoutLeft, Kokkos::HostSpace> mm_u_host(b, dim_host_ptr[0], dim_host_ptr[1]);
            Kokkos::deep_copy(mm_u_host, mm_u);

            // Real* copied_ptr = const_cast<Real*>(mm_u_host.data());
            // printf("--- after deep bis %i ---\n\n", *iter);
        }
        PDI_release("m_u_host_kokkos_view_dimensions");
        PDI_release("m_u_kokkos_view_dimensions");
        PDI_release("m_u_host");
        PDI_release("m_u");
        PDI_release("time");
        PDI_release("freq");
        PDI_release("iter");
    }
}

void GodunovSolver::pdiExposeData()
{
    Kokkos::fence();
    std::chrono::steady_clock::time_point m_start_io = std::chrono::steady_clock::now();

    printf("*********** ici %s ***********************\n", m_params->output.type.c_str());
    if(m_params->output.type.c_str() == "gpu_pdi") 
    {
    printf("*********** if ***********************\n");
 
#if defined(Euler_ENABLE_PDI)
    std::array<size_t, 2> m_u_kokkos_view_dimensions = { m_u.extent(0), m_u.extent(1) };
    std::array<size_t, 2> m_u_host_kokkos_view_dimensions = { m_u_host.extent(0), m_u_host.extent(1) };
 
    std::array<int, 3> pdi_ncells;
    pdi_ncells[IX] = m_grid.m_nbCells[IX] * m_grid.m_dom[IX];
    pdi_ncells[IY] = m_grid.m_nbCells[IY] * m_grid.m_dom[IY];
    pdi_ncells[IZ] = m_grid.m_nbCells[IZ] * m_grid.m_dom[IZ];
 
    std::array<int, 3> pdi_ncells_local;
    pdi_ncells_local[IX] = m_grid.m_nbCells[IX];
    pdi_ncells_local[IY] = m_grid.m_nbCells[IY];
    pdi_ncells_local[IZ] = m_grid.m_nbCells[IZ];
 
    int tmp_rank=0;
#if defined(MPI_SESSION)
    MPI_Comm_rank(MPI_COMM_WORLD, &tmp_rank);
    std::array<int, three_d> mpi_coords = m_grid.comm.getCoords(m_grid.comm.rank());
#endif
 
    std::array<int, 3> pdi_start;
    pdi_start[IX] = m_grid.m_nbCells[IX] * mpi_coords[IX];
    pdi_start[IY] = m_grid.m_nbCells[IY] * mpi_coords[IY];
    pdi_start[IZ] = m_grid.m_nbCells[IZ] * mpi_coords[IZ];

    char *prefix_c_str;
    PDI_access("prefix", (void **)&prefix_c_str, PDI_IN);
    // printf("prefix_c_str %s \n", prefix_c_str);
    std::string prefix(prefix_c_str);
    // printf("prefix %s \n", prefix.c_str());
    PDI_release("prefix");

    std::array<Real, 3> origin;
    origin[IX] = m_grid.m_lowGlobal[IX];
    origin[IY] = m_grid.m_lowGlobal[IY];
    origin[IZ] = m_grid.m_lowGlobal[IZ];
 
    std::array<Real, 3> dl;
    dl[IX] = m_grid.m_dl[IX];
    dl[IY] = m_grid.m_dl[IY];
    dl[IZ] = m_grid.m_dl[IZ];

    PDI_multi_expose("data_GPU_before",
            "iStep", (void*)&(Super::m_iteration), PDI_OUT,
            "m_u", (void*)(m_u.data()), PDI_OUT,
            "m_u_host", (void*)(m_u_host.data()), PDI_OUT,
            "m_u_kokkos_view_dimensions", (void*)&m_u_kokkos_view_dimensions, PDI_OUT,
            "m_u_host_kokkos_view_dimensions", (void*)&m_u_host_kokkos_view_dimensions, PDI_OUT,
            NULL);

    m_writer->write(m_u_host, m_grid, Super::m_iteration, Super::m_t,
                    m_params->thermo.gamma, m_params->thermo.mmw);
#endif

    }
    else
    {
    printf("*********** else ***********************\n");
    #if defined(Euler_ENABLE_PDI)
    PDI_multi_expose("data_on_GPU",
                     "iStep", (void*)&(Super::m_iteration), PDI_OUT,
                     "time", (void*)&(m_t), PDI_OUT,
                     NULL);
    #endif

    if (m_should_save)
    {
       Kokkos::Profiling::pushRegion("I/O - Checkpoint");
       if(Super::m_iteration%100 == 0) Print() << "===================== output at iteration = " << Super::m_iteration << " time t = "<<Super::m_t<< std::endl;
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
