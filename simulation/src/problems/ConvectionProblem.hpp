#pragma once

#include "ConvectionParams.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct ConvectionProblem : Problem
{
    using Array = Array3d;

    ConvectionProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;
    void make_boundaries_user(Array u, const UniformGrid& grid, int idim, int iside) final;

    std::shared_ptr<ConvectionParams> m_prob_params;
};

}}
