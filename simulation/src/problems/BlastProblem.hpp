#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "BlastInitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct BlastProblem : Problem
{
    using Array = Array3d;

    BlastProblem(const std::shared_ptr<Params>& params)
        : Problem {params}
    {
    }

    void initialize(Array u, const UniformGrid& grid) const final
    {
        Print() << "Initializing Blast problem... ";

        BlastInitKernel::apply(*Problem::m_params, u, grid);

        Print() << "done\n";
    }
};


}}
