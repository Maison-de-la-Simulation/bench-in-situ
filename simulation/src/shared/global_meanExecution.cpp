#include "global_meanExecution.hpp"

#include "global_meanKernel.hpp"

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

void
Executeglobal_mean( const Params& params, const UniformGrid& grid, Array3d q, global_mean& sum )
{
    global_meanKernel kernel( params, grid, q );
    global_mean sum_loc(0.);
    Kokkos::MDRangePolicy<Kokkos::Rank<3>> range({2,2,2},{params.mesh.nbCells[IX]+2,params.mesh.nbCells[IY]+2,params.mesh.nbCells[IZ]+2});
    Kokkos::parallel_reduce( "Global mean kernel - MDRangePolicy", range, kernel, sum_loc );

    sum_loc.m_emag /= params.mesh.nbCells[IX]*params.mesh.nbCells[IY]*params.mesh.nbCells[IZ];
    sum_loc.m_ekin /= params.mesh.nbCells[IX]*params.mesh.nbCells[IY]*params.mesh.nbCells[IZ];

#if defined( MPI_SESSION )
    // AllReduce on emag using local copy "emag_loc"
    const int type = grid.comm.template dataType< Real >();

    grid.comm.allReduce( &sum_loc.m_emag, &sum.m_emag, 1, type, grid.comm.SUM );
    grid.comm.allReduce( &sum_loc.m_ekin, &sum.m_ekin, 1, type, grid.comm.SUM );

    sum.m_emag /= grid.comm.size();
    sum.m_ekin /= grid.comm.size();
#else
    sum = sum_loc;
#endif
}

} // namespace hydro
