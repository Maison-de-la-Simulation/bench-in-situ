#include "PrimitiveToConservativeExecution.hpp"

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "PrimitiveToConservativeKernel.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

void
ExecutePrimitiveToConservative( const Params& params, const UniformGrid& grid,
                                const Array3d& q, const Array3d& u )
{
    PrimitiveToConservativeKernel kernel( params, grid, q, u );
    Kokkos::RangePolicy< Int > range( 0, grid.nbCells() );
    Kokkos::parallel_for( "P2C kernel", range, kernel );
}

} // namespace hydro
