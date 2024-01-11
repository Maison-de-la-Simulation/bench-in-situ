#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

void ExecuteConvectionSourceTerm( const Params& params, const UniformGrid& grid,
                         const Array3d& u, const Array3d& q, Real dt );

} // namespace hydro
