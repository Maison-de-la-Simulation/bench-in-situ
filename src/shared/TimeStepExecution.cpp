#include "TimeStepExecution.hpp"

#include "TimeStepKernel.hpp"

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

TimeStep
ExecuteTimeStep( const Params& params, const UniformGrid& grid, Array3d q )
{
    TimeStepKernel kernel( params, grid, q );
    TimeStep invDt;
    Kokkos::RangePolicy< Int > range( 0, grid.nbCells() );
    Kokkos::parallel_reduce( "Time step kernel - RangePolicy", range, kernel, invDt );
    TimeStep dt = TimeStep( std::numeric_limits< Real >::infinity() );
    if ( params.hydro.hydro_enabled )
    {
        dt.hydro = params.run.cfl / invDt.hydro;
    }
#if defined( MPI_SESSION )
    // AllReduce on dt using local copy "dt_loc"
    TimeStep dt_loc = dt;
    const int type = grid.comm.template dataType< Real >();
    if ( params.hydro.hydro_enabled )
    {
        grid.comm.allReduce( &dt_loc.hydro, &dt.hydro, 1, type, grid.comm.MIN );
    }

#endif
    return dt;
}

} // namespace hydro
