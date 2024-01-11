#include "HydroEngine.hpp"

#include "DistributedMemorySession.hpp"
#include "HydroParams.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "ProblemFactory.hpp"
#include "SolverFactory.hpp"

#include <cstring>
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace hydro
{

EngineNd::EngineNd(const std::string& filename)
    : Engine  {}
    , params  {std::make_shared<Params>(filename)}
    , problem {problems::ProblemFactory::New(params->problem, params)}
    , solver  {SolverFactory::New(params->run.solver, problem)}
{
}

void EngineNd::print_configuration(std::ostream& os) const
{
    params->print(os);
    os << std::left  << std::setw(40) << std::setfill('.') << "memory requested [in Mo]"
        << std::right << std::setw(40) << std::setfill('.') << solver->memoryUsage() << '\n';
    os << "Starting simulation from ";
    os << std::setprecision(std::numeric_limits<Real>::digits10);
    os << std::scientific;
    os << "step n=";
    os << std::setw(std::numeric_limits<int>::digits10) << std::setfill('.') << solver->iteration();
    os << "; time t=" << solver->time() << "\n";
}

void EngineNd::run() const
{
    Kokkos::Profiling::pushRegion("I/O");
    solver->pdiExposeData();
    Kokkos::Profiling::popRegion();

    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    Kokkos::Profiling::pushRegion("Time loop");

    bool safe_save_after_20h=true;
    std::chrono::hours time_safe_save(20);
    int freq_check_safe_save=100;

    while (!solver->finished())
    {
        Kokkos::Profiling::pushRegion("Time step");
        Real dt {solver->computeTimeStep()};
        Kokkos::Profiling::popRegion();

        solver->prepareNextOutput(dt);

        Kokkos::Profiling::pushRegion("Hydrodynamical scheme");
        solver->nextIteration(dt);
        Kokkos::Profiling::popRegion();

        if (solver->shouldPrintInformation())
        {
            solver->printInformation(dt);
        }

        Kokkos::Profiling::pushRegion("I/O");

        std::chrono::steady_clock::time_point date = std::chrono::steady_clock::now();
        std::chrono::steady_clock::duration duration = date-start;

        if (solver->iteration()%freq_check_safe_save==0)
        {
            bool should_save_and_quit=(safe_save_after_20h)&&(duration>=time_safe_save);

#if defined(MPI_SESSION)
            MPI_Bcast(&should_save_and_quit, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif

            if (should_save_and_quit)
            {
              std::cout<<"20h of simulation reached, saving solution and stopping simulation ..."<<std::endl;
              solver->set_should_save();
              solver->set_time_limit_reached();
            }
        }
        solver-> pdiExposeData();
        Kokkos::Profiling::popRegion();
    }

    Kokkos::Profiling::popRegion();
    const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    solver->printMonitoring(std::chrono::duration<double>(end-start).count());
}

std::shared_ptr<Engine> EngineFactory::New(const std::string& filename)
{
    const INIReader reader {filename};
    if(reader.ParseError()<0)
    {
        throw std::runtime_error("Error opening file \""+filename+"\": "+std::strerror(errno));
    }

    std::shared_ptr<Engine> engine = nullptr;
    int dim {reader.GetInteger("problem", "dimension", 0)};
    if (dim == three_d)
    {
        engine = std::make_shared<EngineNd>(filename);
    }
    else
    {
        throw std::runtime_error("Wrong dimension. Dimension allowed: 3");
    }

    return engine;
}

}  // hydro
