#pragma once

#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "RayleighTaylorParams.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct RayleighTaylorProblem : Problem
{
    using Array = Array3d;

    RayleighTaylorProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;
    void make_boundaries_user(Array u, const UniformGrid& grid, int idim, int iside) final;

    std::shared_ptr<RayleighTaylorParams> m_problem_params;
};

}}
