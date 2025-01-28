#include "MusclReconstructionExecution.hpp"
#include "MusclReconstructionKernel.hpp"

namespace hydro {

void ExecuteMusclReconstruction(const Params &params, const UniformGrid &grid,
                                Array3d q,
                                const Kokkos::Array<Array3d, 2 * three_d> &qr,
                                Real dt)
{
    using Kernel = MusclReconstructionKernel;
    Kernel kernel(params, grid, q, qr, dt);
    Kokkos::parallel_for("Muscl reconstruction kernel - MDRangePolicy",
                         Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                         {2 - Kernel::ghostDepth, 2 - Kernel::ghostDepth, 2 - Kernel::ghostDepth},
                         {grid.m_nbCells[ IX ] + 2 + Kernel::ghostDepth, grid.m_nbCells[ IY ] + 2 + Kernel::ghostDepth, grid.m_nbCells[ IZ ] + 2 + Kernel::ghostDepth},
                         {64, 2, 2}),
                         kernel);
}

} // namespace hydro
