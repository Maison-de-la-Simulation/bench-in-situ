#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroTypes.hpp"
#include "HydroParams.hpp"

#include <cmath>
#include <Kokkos_Core.hpp>

namespace hydro
{

void ExecutePrimitiveToConservative( const Params& params, const UniformGrid& grid,
                                     const Array3d& q, const Array3d& u );

}  // hydro
