#include "OTangProblem.hpp"

#include "OTangInitKernel.hpp"
#include "OTangParams.hpp"
#include "HydroParams.hpp"

#include <memory>

namespace hydro { namespace problems
{

OTangProblem::OTangProblem(const std::shared_ptr<Params>& params)
    : Problem   {params}
    , m_problem_params {std::make_shared<OTangParams>(params->reader)}
{
}

void OTangProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing OTang problem... ";

    OTangInitKernel::apply(*m_params, *m_problem_params, u, grid);

    Print() << "done\n";
}

}}
