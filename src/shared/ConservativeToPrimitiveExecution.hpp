#pragma once

#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

namespace hydro
{

void ExecuteConservativeToPrimitive( const Params& params, const UniformGrid& grid,
                                     const Array3d& u, const Array3d& q );

} // namespace hydro
