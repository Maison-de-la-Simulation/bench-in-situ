#pragma once

#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

void FluxesAndUpdateKernelDispatch( const Params& params, const UniformGrid& grid,
                                    const Array3d& u,
                                    const Array3d& q,
                                    const Kokkos::Array<Array3d, 2*3>& qr,
                                    Real dt );

} // namespace hydro
