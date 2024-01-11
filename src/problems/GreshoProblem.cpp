#include "GreshoProblem.hpp"

#include "GreshoInitKernel.hpp"
#include "GreshoParams.hpp"
#include "HydroParams.hpp"

#include <memory>

namespace hydro { namespace problems
{

GreshoProblem::GreshoProblem(const std::shared_ptr<Params>& params)
    : Problem   {params}
    , m_problem_params {std::make_shared<GreshoParams>(params->reader)}
{
}

void GreshoProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing Gresho problem... ";

    GreshoInitKernel::apply(*m_params, *m_problem_params, u, grid);

    Print() << "done\n";
}

}}
