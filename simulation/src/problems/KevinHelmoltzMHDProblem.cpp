#include "HydroProblem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "KevinHelmoltzMHDProblem.hpp"
#include "KevinHelmoltzMHDInitKernel.hpp"
#include "KevinHelmoltzMHDBoundariesKernel.hpp"


#include <memory>

namespace hydro { namespace problems
{

KevinHelmoltzMHDProblem::KevinHelmoltzMHDProblem(const std::shared_ptr<Params>& params)
    : Problem {params}
{
}

void KevinHelmoltzMHDProblem::initialize(Array u, const UniformGrid& grid) const
{
    Print() << "Initializing KevinHelmoltzMHD problem... ";

    KevinHelmoltzMHDInitKernel::apply(*m_params, u, grid);

    Print() << "done\n";
}

void KevinHelmoltzMHDProblem::make_boundaries_user(Array u, const UniformGrid& grid, int idim, int iside)
{
    KevinHelmoltzMHDBoundariesKernel::apply(u, grid, *m_params, idim, iside);
}


}}
