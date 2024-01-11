#pragma once

#include "OTangParams.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"
#include "inih/INIReader.hpp"

#include <memory>

namespace hydro { namespace problems
{

struct OTangProblem : Problem
{
    using Array = Array3d;

    OTangProblem(const std::shared_ptr<Params>& params);
    void initialize(Array u, const UniformGrid& grid) const final;

    std::shared_ptr<OTangParams> m_problem_params;
};

}}
