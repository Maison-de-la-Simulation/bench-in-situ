#include "FluxesAndUpdateKernelDispatch.hpp"

#include "FluxesAndUpdateExecution.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "ars/HydroRiemannAllRegime.hpp"
#include "ars/MHDRiemann_AllRegime_3W.hpp"
#include "ars/MHDRiemann_AllRegime_5W.hpp"
#include "ars/MHDRiemann_ARWB_5W.hpp"
#include "ars/MHD3W_optimized.hpp"
#include "ars/HydroRiemannAllRegime_experiment_entropy.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

void FluxesAndUpdateKernelDispatch( const Params& params, const UniformGrid& grid,
                                    const Array3d& u,
                                    const Array3d& q,
                                    const Kokkos::Array<Array3d, 2*3>& qr,
                                    Real dt )
{
   if ( params.run.riemann == "all_regime" )
    {
        ExecuteFluxesAndUpdate< riemann::AllRegime >( params, grid, u, q, qr,  dt );
    }
    else if ( params.run.riemann == "mhd_all_regime_3w" )
    {
        ExecuteFluxesAndUpdate< riemann::MHD_AllRegime_3W>( params, grid, u, q, qr,  dt );
    }
    else if ( params.run.riemann == "mhd_all_regime_5w" )
    {
        ExecuteFluxesAndUpdate< riemann::MHD_AllRegime_5W>( params, grid, u, q, qr,  dt );
    }
    else if ( params.run.riemann == "mhd_ARWB_5w" )
    {
        ExecuteFluxesAndUpdate< riemann::MHD_ARWB_5W>( params, grid, u, q, qr,  dt );
    }
    else if ( params.run.riemann == "MHD3W_optimized" )
    {
        ExecuteFluxesAndUpdate< riemann::MHD3W_optimized>( params, grid, u, q, qr,  dt );
    }
    else if ( params.run.riemann == "exp_entropy" )
    {
        ExecuteFluxesAndUpdate< riemann::AllRegime_exp>( params, grid, u, q, qr,  dt );
    }
    else
    {
        throw std::runtime_error("Unknown Riemann solver");
    }
}

} // namespace hydro
