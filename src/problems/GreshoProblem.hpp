#pragma once

#include "GreshoParams.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "inih/INIReader.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct GreshoProblem : Problem
{
    using Array = Array3d;

    GreshoProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;

    std::shared_ptr<GreshoParams> m_problem_params;
};

}}
