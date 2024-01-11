#pragma once

#include "HydroSolver.hpp"
#include "HydroProblem.hpp"

#include <memory>
#include <string>

namespace hydro
{

class SolverFactory
{
public:
    static std::shared_ptr<Solver> New(const std::string& solver_name,
                                       std::shared_ptr<Problem> problem);
};


}
