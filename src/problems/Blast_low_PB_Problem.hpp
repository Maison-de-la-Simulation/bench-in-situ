#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "Blast_low_PB_InitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct Blast_low_PB_Problem : Problem
{
    using Array = Array3d;

    Blast_low_PB_Problem(const std::shared_ptr<Params>& params)
        : Problem {params}
    {
    }

    void initialize(Array u, const UniformGrid& grid) const final
    {
        Print() << "Initializing Blast_low_PB_ problem... ";

        Blast_low_PB_InitKernel::apply(*Problem::m_params, u, grid);

        Print() << "done\n";
    }
};


}}
