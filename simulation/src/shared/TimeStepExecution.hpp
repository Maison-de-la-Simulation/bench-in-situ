#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "TimeStep.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

TimeStep ExecuteTimeStep( const Params& params, const UniformGrid& grid,
                          Array3d q );
                          
} // namespace hydro
