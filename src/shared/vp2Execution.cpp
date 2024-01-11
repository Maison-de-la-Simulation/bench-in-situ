#include "vp2Execution.hpp"

#include "slice_averageKernel.hpp"

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "slice_average.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace hydro
{

void
Executevp2(const Params& params, const UniformGrid& grid, Array3d q, Kokkos::View<double**,Kokkos::HostSpace> profiles)
{
    const int nx = params.mesh.nbCells[IX];
    const int ny = params.mesh.nbCells[IY];
    const int nz = params.mesh.nbCells[IZ];

    for (Int iz = 2; iz < nz+2; iz++)
    {
        slice_averageKernel kernel( params, grid, q, iz );
        slice_average slice(0.0);
        Kokkos::MDRangePolicy<Kokkos::Rank<2>> range({2,2},{nx+2,ny+2});
        Kokkos::parallel_reduce( "Slice average kernel - MDRangePolicy", range, kernel, slice );

        profiles(0,iz-2) = slice.m_Bx2/(nx*ny);
        profiles(1,iz-2) = slice.m_By2/(nx*ny);
        profiles(2,iz-2) = slice.m_Bz2/(nx*ny);
        profiles(3,iz-2) = slice.m_u/(nx*ny);
        profiles(4,iz-2) = slice.m_v/(nx*ny);
        profiles(5,iz-2) = slice.m_logmu/(nx*ny);
        profiles(6,iz-2) = slice.m_logtheta/(nx*ny);
    }
}

} // namespace hydro
