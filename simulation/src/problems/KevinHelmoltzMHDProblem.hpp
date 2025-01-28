#pragma once

#include "HydroProblem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct KevinHelmoltzMHDProblem : Problem
{
    using Array = Array3d;

    KevinHelmoltzMHDProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;
    void make_boundaries_user(Array u, const UniformGrid& grid, int idim, int iside) final;

};

}}
