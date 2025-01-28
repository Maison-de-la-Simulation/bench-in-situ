#pragma once

#include "HydroProblem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct SkewedShockProblem : Problem
{
    using Array = Array3d;

    SkewedShockProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;
};

}}
