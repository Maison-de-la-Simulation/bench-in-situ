#include "ConservativeToPrimitiveExecution.hpp"

#include "ConservativeToPrimitiveKernel.hpp"

namespace hydro
{

void
ExecuteConservativeToPrimitive( const Params& params, const UniformGrid& grid,
                                const Array3d& u, const Array3d& q )
{
    ConservativeToPrimitiveKernel kernel( params, grid, u, q );
    Kokkos::RangePolicy< Int > range( 0, grid.nbCells() );
    Kokkos::parallel_for( "C2P kernel", range, kernel );
}

} // namespace hydro
