#include "ConvectionProblem.hpp"

#include "ConvectionBoundariesKernel.hpp"
#include "ConvectionInitKernel.hpp"
#include "ConvectionParams.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroTypes.hpp"

#include <memory>

namespace hydro { namespace problems
{

ConvectionProblem::ConvectionProblem(const std::shared_ptr<Params>& params)
    : Problem {params}
    , m_prob_params  {std::make_shared<ConvectionParams>(params->reader)}
{
}

void ConvectionProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing Convection problem... ";

    ConvectionInitKernel::apply(*m_params, *m_prob_params, u, grid);

    Print() << "done" << std::endl;
}

void ConvectionProblem::make_boundaries_user(Array u, const UniformGrid& grid, int, int iside)
{
    ConvectionBoundariesKernel::apply(u, grid, *m_params, *m_prob_params, iside);
}

}}
