#pragma once

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
void ExecuteFluxesAndUpdate( const Params& params, const UniformGrid& grid,
                             const Array3d& u, const Array3d& q, const Kokkos::Array<Array3d, 2*three_d>& qr,Real dt );

#define ETI_ExecuteFluxesAndUpdate( ars )                                                     \
    void ExecuteFluxesAndUpdate< ars >(                                                       \
        const Params& params, const UniformGrid& grid, const Array3d& u,      \
        const Array3d& q, const Kokkos::Array<Array3d, 2*three_d>& qr, Real dt )


extern template ETI_ExecuteFluxesAndUpdate( riemann::AllRegime );

extern template ETI_ExecuteFluxesAndUpdate( riemann::MHD_AllRegime_3W );

extern template ETI_ExecuteFluxesAndUpdate( riemann::MHD_AllRegime_5W );

extern template ETI_ExecuteFluxesAndUpdate( riemann::MHD_ARWB_5W );

extern template ETI_ExecuteFluxesAndUpdate( riemann::MHD3W_optimized );

extern template ETI_ExecuteFluxesAndUpdate( riemann::AllRegime_exp );


} // namespace hydro
