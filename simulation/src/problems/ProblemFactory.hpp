#pragma once

#include "HydroTypes.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"

#include <memory>

namespace hydro { namespace problems
{

// Be careful, here we do not use explicit template instantiation but
// full specialization.
class ProblemFactory
{
public:
    static std::shared_ptr<Problem> New(const std::string& problem_name,
                                             const std::shared_ptr<Params>& params);
};

}}
