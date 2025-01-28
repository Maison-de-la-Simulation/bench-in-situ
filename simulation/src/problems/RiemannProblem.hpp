#pragma once

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "RiemannInitKernel.hpp"
#include "RiemannParams.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct RiemannProblem : Problem
{
    using Array = Array3d;

    RiemannProblem(const std::shared_ptr<Params>& params)
        : Problem   {params}
        , m_problem_params {std::make_shared<RiemannParams>(params->reader)}
    {
    }

    void initialize(Array u, const UniformGrid& grid) const final
    {
        Print() << "Initializing shock tube problem... ";

        RiemannInitKernel::apply(*Problem::m_params, *m_problem_params, u, grid);

        Print() << "done\n";
    }

    std::shared_ptr<RiemannParams> m_problem_params;
};

}}
