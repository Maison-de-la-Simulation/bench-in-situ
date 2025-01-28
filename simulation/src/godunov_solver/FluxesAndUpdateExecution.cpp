#include "FluxesAndUpdateExecution.hpp"
#include "GodunovFluxesKernel.hpp"

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
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

template < class RiemannSolver >
void
ExecuteFluxesAndUpdate( const Params& params, const UniformGrid& grid,
                        const Array3d& u, const Array3d& q,  const Kokkos::Array<Array3d, 2*three_d>& qr, Real dt )
{
    using Kernel = FluxesAndUpdateKernel< RiemannSolver >;
    Kernel kernel( params, grid, u, q, qr, dt );
    Kokkos::parallel_for("Godunov kernel - MDRangePolicy",
                         Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                         {2 - Kernel::ghostDepth, 2 - Kernel::ghostDepth, 2 - Kernel::ghostDepth},
                         {grid.m_nbCells[ IX ] + 2 + Kernel::ghostDepth, grid.m_nbCells[ IY ] + 2 + Kernel::ghostDepth, grid.m_nbCells[ IZ ] + 2 + Kernel::ghostDepth},
                         {64, 2, 2}),
                         kernel);

}

template ETI_ExecuteFluxesAndUpdate(riemann::AllRegime );

template ETI_ExecuteFluxesAndUpdate( riemann::MHD_AllRegime_3W );

template ETI_ExecuteFluxesAndUpdate( riemann::MHD_AllRegime_5W );

template ETI_ExecuteFluxesAndUpdate( riemann::MHD_ARWB_5W );

template ETI_ExecuteFluxesAndUpdate( riemann::MHD3W_optimized );

template ETI_ExecuteFluxesAndUpdate(riemann::AllRegime_exp );

} // namespace hydro
