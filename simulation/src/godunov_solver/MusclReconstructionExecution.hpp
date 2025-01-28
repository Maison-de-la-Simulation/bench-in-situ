#pragma once

#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

void ExecuteMusclReconstruction(const Params &params, const UniformGrid &grid,
                                Array3d q,
                                const Kokkos::Array<Array3d, 2 * three_d> &qr,
                                Real dt);

} // namespace hydro
