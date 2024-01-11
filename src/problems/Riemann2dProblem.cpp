#include "HydroProblem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "Riemann2dProblem.hpp"
#include "Riemann2dInitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

Riemann2dProblem::Riemann2dProblem(const std::shared_ptr<Params>& params)
    : Problem {params}
{
}

void Riemann2dProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing Riemann2d problem... ";

    Riemann2dInitKernel::apply(*m_params, u, grid);

    Print() << "done\n";
}

}}
