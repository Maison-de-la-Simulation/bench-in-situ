#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "ImplodeInitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct ImplodeProblem : Problem
{
    using Array = Array3d;

    ImplodeProblem(const std::shared_ptr<Params>& params)
        : Problem {params}
    {
    }

    void initialize(Array u, const UniformGrid& grid) const final
    {
        Print() << "Initializing implode problem... ";

        ImplodeInitKernel::apply(*Problem::m_params, u, grid);

        Print() << "done\n";
    }
};


}}
