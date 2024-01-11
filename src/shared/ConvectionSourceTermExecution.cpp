#include "ConvectionSourceTermExecution.hpp"

#include "ConvectionSourceTermKernel.hpp"
#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

namespace hydro
{

void
ExecuteConvectionSourceTerm( const Params& params, const UniformGrid& grid,
                    const Array3d& u, const Array3d& q, Real dt )
{
    using Kernel = ConvectionSourceTermKernel;
    Kernel kernel(params, grid, u, q, dt );
    Kokkos::parallel_for("Convection source term kernel - MDRangePolicy",
                         Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                         {2 - Kernel::ghostDepth, 2 - Kernel::ghostDepth, 2 - Kernel::ghostDepth},
                         {grid.m_nbCells[ IX ] + 2 + Kernel::ghostDepth, grid.m_nbCells[ IY ] + 2 + Kernel::ghostDepth, grid.m_nbCells[ IZ ] + 2 + Kernel::ghostDepth},
                         {64, 2, 2}),
                         kernel);
}

} // namespace hydro
