#include "SolverFactory.hpp"
#include "godunov_solver/Godunov.hpp"
#include "HydroSolver.hpp"
#include "HydroProblem.hpp"

#include <exception>
#include <memory>
#include <string>

namespace hydro
{

std::shared_ptr<Solver> SolverFactory::New(const std::string& solver_name,
                                                std::shared_ptr<Problem> problem)
{
    std::shared_ptr<Solver> solver = nullptr;
    if (solver_name == "godunov")
    {
        solver = std::make_shared<GodunovSolver>(problem);
    }
    else
    {
        throw std::runtime_error("Unknown solver");
    }

    return solver;
}

}
