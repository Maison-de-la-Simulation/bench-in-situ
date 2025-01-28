#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

void Executevp2( const Params& params, const UniformGrid& grid,
                          Array3d q, Kokkos::View<double**,Kokkos::HostSpace> profiles );

} // namespace hydro
