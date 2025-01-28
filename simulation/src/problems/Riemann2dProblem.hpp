#pragma once

#include "HydroProblem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct Riemann2dProblem : Problem
{
    using Array = Array3d;

    Riemann2dProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;
};

}}
