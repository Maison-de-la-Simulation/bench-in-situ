#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "global_mean.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

void Executeglobal_mean( const Params& params, const UniformGrid& grid,
                         Array3d q, global_mean& sum );

} // namespace hydro
