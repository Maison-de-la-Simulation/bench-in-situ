#include "HydroProblem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "SkewedShockProblem.hpp"
#include "SkewedShockInitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

SkewedShockProblem::SkewedShockProblem(const std::shared_ptr<Params>& params)
    : Problem {params}
{
}

void SkewedShockProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing SkewedShock problem... ";

    SkewedShockInitKernel::apply(*m_params, u, grid);

    Print() << "done\n";
}

}}
