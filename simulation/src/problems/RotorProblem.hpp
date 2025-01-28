#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "RotorInitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct RotorProblem : Problem
{
    using Array = Array3d;

    RotorProblem(const std::shared_ptr<Params>& params)
        : Problem {params}
    {
    }

    void initialize(Array u, const UniformGrid& grid) const final
    {
        Print() << "Initializing Rotor problem... ";

        RotorInitKernel::apply(*Problem::m_params, u, grid);

        Print() << "done\n";
    }
};


}}
