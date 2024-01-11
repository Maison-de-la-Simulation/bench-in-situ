#include "ReaderPDI.hpp"

#include "DistributedMemorySession.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "ReaderBase.hpp"
#include "MHDSystem.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <pdi.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <assert.h>

namespace hydro { namespace io
{

ReaderPDI::ReaderPDI(const UniformGrid&, const Params& params,
                     const std::vector<std::pair<int, std::string>>& variables)
    : filename {params.run.restart_filename}
    , m_variables {variables}
{
}


void ReaderPDI::read(HostArrayDyn u, const UniformGrid& grid,
                             Int& iStep, Real& time, Int& outputId, Int& restartId)
{
    // inquery the restart grid size
    int grid_size_read[3];
    int filename_size = filename.size();

    PDI_multi_expose("prep",
                     "grid_size", &grid_size_read, PDI_IN,
                     "restart_filename_size", &filename_size, PDI_INOUT,
                     "restart_filename", filename.data(), PDI_INOUT,
                     NULL);

    assert((grid.m_nbCells[IX] * grid.m_dom[IX]) % grid_size_read[IX] == 0);
    assert((grid.m_nbCells[IY] * grid.m_dom[IY]) % grid_size_read[IY] == 0);
    assert((grid.m_nbCells[IZ] * grid.m_dom[IZ]) % grid_size_read[IZ] == 0);

    int ratio_interpolation[3];
    ratio_interpolation[IX] = grid.m_nbCells[IX] * grid.m_dom[IX] / grid_size_read[IX];
    ratio_interpolation[IY] = grid.m_nbCells[IY] * grid.m_dom[IY] / grid_size_read[IY];
    ratio_interpolation[IZ] = grid.m_nbCells[IZ] * grid.m_dom[IZ] / grid_size_read[IZ];

    int ratio=ratio_interpolation[IX]*ratio_interpolation[IY]*ratio_interpolation[IZ];

    PDI_multi_expose("",
                     "ratioX", &ratio_interpolation[IX], PDI_INOUT,
                     "ratioY", &ratio_interpolation[IY], PDI_INOUT,
                     "ratioZ", &ratio_interpolation[IZ], PDI_INOUT,
                     NULL);

    if( ratio == 1) // if same, proceed as usual
    {
        std::cout<<"==== no interpolation needed ===="<<std::endl;
        PDI_multi_expose("read",
                     "iStep", &iStep, PDI_INOUT,
                     "grid_size", &grid_size_read, PDI_IN,
                     "time", &time, PDI_INOUT,
                     "output_id", &outputId, PDI_IN,
                     "restart_id", &restartId, PDI_IN,
                     "local_full_field", u.data(), PDI_INOUT,
                     NULL);
    }

    else  //if not, read on a low resolution grid tmp_field,
          //then interpolate on high resolution grid
    {
        std::cout<<"==== perform spatial interpolation ===="<<std::endl;
        #if defined(MPI_SESSION)
        std::array<int, three_d> coords {grid.comm.getCoords(grid.comm.rank())};
        #else
        std::array<int, three_d> coords {};
        #endif

        std::array<int, 3> read_start;
        read_start[IX] = grid_size_read[IX]/grid.m_dom[IX] * coords[IX];
        read_start[IY] = grid_size_read[IY]/grid.m_dom[IY] * coords[IY];
        read_start[IZ] = grid_size_read[IZ]/grid.m_dom[IZ] * coords[IZ];

        PDI_multi_expose("init_read", "read_start", read_start.data(), PDI_INOUT, NULL);

        Kokkos::View<double****, Layout, Kokkos::HostSpace> field_read("",
                                                                       grid_size_read[IX]/grid.m_dom[IX],
                                                                       grid_size_read[IY]/grid.m_dom[IY],
                                                                       grid_size_read[IZ]/grid.m_dom[IZ],
                                                                       MHDSystem::nbvar);

        PDI_multi_expose("event_read_fields",
                     "iStep", &iStep, PDI_INOUT,
                     "time", &time, PDI_INOUT,
                     "output_id", &outputId, PDI_IN,
                     "restart_id", &restartId, PDI_IN,
                     "d", &field_read(0, 0, 0, MHDSystem::VarCons::ID), PDI_INOUT,
                     "E", &field_read(0, 0, 0, MHDSystem::VarCons::IE), PDI_INOUT,
                     "mx", &field_read(0, 0, 0, MHDSystem::VarCons::IMx), PDI_INOUT,
                     "my", &field_read(0, 0, 0, MHDSystem::VarCons::IMy), PDI_INOUT,
                     "mz", &field_read(0, 0, 0, MHDSystem::VarCons::IMz), PDI_INOUT,
                     "Bx", &field_read(0, 0, 0, MHDSystem::VarCons::IBx), PDI_INOUT,
                     "By", &field_read(0, 0, 0, MHDSystem::VarCons::IBy), PDI_INOUT,
                     "Bz", &field_read(0, 0, 0, MHDSystem::VarCons::IBz), PDI_INOUT,
                     "dX", &field_read(0, 0, 0, MHDSystem::VarCons::I_dX), PDI_INOUT,
                     NULL);

        // interpolate
        int const nx_glo = grid.m_nbCells[IX]+2*grid.m_ghostWidths[IX];
        int const ny_glo = grid.m_nbCells[IY]+2*grid.m_ghostWidths[IY];
        int const nz_glo = grid.m_nbCells[IZ]+2*grid.m_ghostWidths[IZ];

        Kokkos::View<Real ****, Layout, Kokkos::HostSpace> u4d(u.data(), nx_glo, ny_glo, nz_glo, MHDSystem::nbvar);

        for(int ivar=0; ivar<MHDSystem::nbvar; ivar++)
        {
            Kokkos::parallel_for("", Kokkos::MDRangePolicy<Kokkos::Rank<3>, Kokkos::DefaultHostExecutionSpace>({0, 0, 0}, {grid_size_read[IX]/grid.m_dom[IX], grid_size_read[IY]/grid.m_dom[IY], grid_size_read[IZ]/grid.m_dom[IZ]}),
                                 KOKKOS_LAMBDA(int const i, int const j, int const k) {
                                     double const var = field_read(i, j, k, ivar);

                                     int const I = grid.m_ghostWidths[IX]+ratio_interpolation[IX]*i;
                                     int const J = grid.m_ghostWidths[IY]+ratio_interpolation[IY]*j;
                                     int const K = grid.m_ghostWidths[IZ]+ratio_interpolation[IZ]*k;

                                     for ( int ii = 0; ii < ratio_interpolation[IX]; ++ii )
                                     {
                                         for ( int jj = 0; jj < ratio_interpolation[IY]; ++jj )
                                         {
                                             for ( int kk = 0; kk < ratio_interpolation[IZ]; ++kk )
                                             {
                                                 u4d(I+ii, J+jj, K+kk, ivar) = var;
                                             }
                                         }
                                     }
                                 });
        }
    }
}

}}
