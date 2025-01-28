#include "RayleighTaylorProblem.hpp"

#include "RayleighTaylorBoundariesKernel.hpp"
#include "RayleighTaylorInitKernel.hpp"

#include <memory>

namespace hydro { namespace problems
{

RayleighTaylorProblem::RayleighTaylorProblem(const std::shared_ptr<Params>& params)
    : Problem   {params}
    , m_problem_params {std::make_shared<RayleighTaylorParams>(params->reader)}
{
}

void RayleighTaylorProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing Rayleigh-Taylor problem... ";

    RayleighTaylorInitKernel::apply(*m_params, *m_problem_params, u, grid);

    Print() << "done\n";
}

void RayleighTaylorProblem::make_boundaries_user(Array u, const UniformGrid& grid, int idim, int iside)
{
    RayleighTaylorBoundariesKernel::apply(u, grid, *m_params, idim, iside);
}

}}
