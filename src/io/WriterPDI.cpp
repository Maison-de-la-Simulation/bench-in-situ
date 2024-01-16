#include "WriterPDI.hpp"

#include "DistributedMemorySession.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "WriterBase.hpp"
#include "Utils.hpp"
#include "MHDSystem.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <pdi.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace hydro { namespace io
{

std::string getFilename(std::string const &prefix, Int outputId) {
  // write outputId in string outputNum
  std::ostringstream outputNum;
  outputNum << std::setw(std::numeric_limits<Int>::digits10);
  outputNum << std::setfill('0');
  outputNum << outputId;

  // concatenate file prefix + file number + suffix
  std::string filename(prefix);
  filename += "_" + outputNum.str();
  filename += ".h5";
  return filename;
}

extern "C"
{

    void writeXML() {
        if (!Session::isIOProc())
        {
            return;
        }

        Int* restartId;
        PDI_access("restart_id", (void**)&restartId, PDI_IN);
        std::ostringstream restartNum;
        restartNum << std::setw(std::numeric_limits<Int>::digits10);
        restartNum << std::setfill('0');
        restartNum << *restartId;
        PDI_release("restart_id");

        char* prefix_c_str;
        PDI_access("prefix", (void**)&prefix_c_str, PDI_IN);
        std::string prefix(prefix_c_str);
        PDI_release("prefix");
        const std::string xdmfFilenameFull{// directory + '/' +
            prefix + '_' +
            restartNum.str() + ".xmf"};
        std::ofstream xdmfFile(xdmfFilenameFull, std::ofstream::trunc);

        xdmfFile << "<?xml version=\"1.0\"?>\n";
        xdmfFile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        xdmfFile << "<Xdmf Version=\"2.0\">\n";
        xdmfFile << std::string(2, ' ') << "<Domain>\n";
        xdmfFile << std::string(4, ' ');
        xdmfFile << "<Grid";
        xdmfFile << " Name=" << '"' << "TimeSeries" << '"';
        xdmfFile << " GridType=" << '"' << "Collection" << '"';
        xdmfFile << " CollectionType=" << '"' << "Temporal" << '"';
        xdmfFile << ">\n";

        std::pair<Int, Real>* outputs_record;
        int* outputs_record_size;
        PDI_access("outputs_record", (void**)&outputs_record, PDI_IN);
        PDI_access("outputs_record_size", (void**)&outputs_record_size, PDI_IN);

        int *ncells;
        PDI_access("ncell", (void **)&ncells, PDI_IN);

        Real *origin;
        PDI_access("origin", (void **)&origin, PDI_IN);

        Real *dl;
        PDI_access("dl", (void **)&dl, PDI_IN);

        int precision = sizeof(Real);

        std::vector<std::string> var_names = hydro::MHDSystem::cons_names();
        std::vector<std::pair<int, std::string>> variables_to_save;
        for (int ivar = 0; ivar < hydro::MHDSystem::nbvar; ++ivar) {
          variables_to_save.push_back(std::make_pair(ivar, var_names[ivar]));
        }

        for ( std::pair<Int, Real>* it = outputs_record; it != outputs_record+*outputs_record_size; ++it ) {
            xdmfFile << std::string(6, ' ');
            xdmfFile << "<Grid Name=" << '"' << "output" << '"';
            xdmfFile << " GridType=" << '"' << "Uniform" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(8, ' ');
            xdmfFile << "<Time Value=" << '"' << it->second << '"'
                     << "/>\n";

            // topology CoRectMesh
            xdmfFile << std::string(8, ' ');
            xdmfFile << "<Topology";
            xdmfFile << " TopologyType=" << '"' << "3DCoRectMesh" << '"';
            xdmfFile << " Dimensions=" << '"';
            for (int idim = 2; idim >= 0; --idim) {
                xdmfFile << ncells[idim] + 1;
                xdmfFile << (idim == 0 ? "\"" : " ");
            }
            xdmfFile << "/>\n";

            // geometry
            xdmfFile << std::string(8, ' ');
            xdmfFile << "<Geometry";
            xdmfFile << " GeometryType=" << '"' << "ORIGIN_DXDYDZ" << '"';
            xdmfFile << ">\n";

            xdmfFile << std::string(10, ' ');
            xdmfFile << "<DataItem";
            xdmfFile << " Name=" << '"' << "Origin" << '"';
            xdmfFile << " NumberType=" << '"' << "Float" << '"';
            xdmfFile << " Precision=" << '"' << precision << '"';
            xdmfFile << " Dimensions=" << '"' << 3 << '"';
            xdmfFile << " Format=" << '"' << "XML" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(12, ' ');
            for (int idim = 2; idim >= 0; --idim) {
                xdmfFile << origin[idim];
                xdmfFile << (idim == 0 ? "\n" : " ");
            }
            xdmfFile << std::string(10, ' ') << "</DataItem>\n";

            xdmfFile << std::string(10, ' ');
            xdmfFile << "<DataItem";
            xdmfFile << " Name=" << '"' << "Spacing" << '"';
            xdmfFile << " NumberType=" << '"' << "Float" << '"';
            xdmfFile << " Precision=" << '"' << precision << '"';
            xdmfFile << " Dimensions=" << '"' << 3 << '"';
            xdmfFile << " Format=" << '"' << "XML" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(12, ' ');
            for (int idim = 2; idim >= 0; --idim) {
                xdmfFile << dl[idim];
                xdmfFile << (idim == 0 ? "\n" : " ");
            }
            xdmfFile << std::string(10, ' ') << "</DataItem>\n";

            xdmfFile << std::string(8, ' ') << "</Geometry>\n";

            // Write gamma
            xdmfFile << std::string(8, ' ');
            xdmfFile << "<Attribute";
            xdmfFile << " Center=" << '"' << "Grid" << '"';
            xdmfFile << " Name=" << '"' << "gamma" << '"';
            xdmfFile << " AttributeType=" << '"' << "Scalar" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(10, ' ');
            xdmfFile << "<DataItem";
            xdmfFile << " NumberType=" << '"' << "Float" << '"';
            xdmfFile << " Precision=" << '"' << precision << '"';
            xdmfFile << " Dimensions=" << '"' << 1 << '"';
            xdmfFile << " Format=" << '"' << "HDF" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(12, ' ') << getFilename(prefix, it->first)
                     << ":/"
                     << "gamma"
                     << "\n";
            xdmfFile << std::string(10, ' ') << "</DataItem>\n";
            xdmfFile << std::string(8, ' ') << "</Attribute>\n";
            // Write mmw
            xdmfFile << std::string(8, ' ');
            xdmfFile << "<Attribute";
            xdmfFile << " Center=" << '"' << "Grid" << '"';
            xdmfFile << " Name=" << '"' << "mmw" << '"';
            xdmfFile << " AttributeType=" << '"' << "Scalar" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(10, ' ');
            xdmfFile << "<DataItem";
            xdmfFile << " NumberType=" << '"' << "Float" << '"';
            xdmfFile << " Precision=" << '"' << precision << '"';
            xdmfFile << " Dimensions=" << '"' << 1 << '"';
            xdmfFile << " Format=" << '"' << "HDF" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(12, ' ') << getFilename(prefix, it->first)
                     << ":/"
                     << "mmw"
                     << "\n";
            xdmfFile << std::string(10, ' ') << "</DataItem>\n";
            xdmfFile << std::string(8, ' ') << "</Attribute>\n";
            // Write Rstar_h
            xdmfFile << std::string(8, ' ');
            xdmfFile << "<Attribute";
            xdmfFile << " Center=" << '"' << "Grid" << '"';
            xdmfFile << " Name=" << '"' << "Rstar_h" << '"';
            xdmfFile << " AttributeType=" << '"' << "Scalar" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(10, ' ');
            xdmfFile << "<DataItem";
            xdmfFile << " NumberType=" << '"' << "Float" << '"';
            xdmfFile << " Precision=" << '"' << precision << '"';
            xdmfFile << " Dimensions=" << '"' << 1 << '"';
            xdmfFile << " Format=" << '"' << "HDF" << '"';
            xdmfFile << ">\n";
            xdmfFile << std::string(12, ' ') << getFilename(prefix, it->first)
                     << ":/"
                     << "Rstar_h"
                     << "\n";
            xdmfFile << std::string(10, ' ') << "</DataItem>\n";
            xdmfFile << std::string(8, ' ') << "</Attribute>\n";

            for (const auto &var : variables_to_save) {
                const std::string var_name = var.second;

                xdmfFile << std::string(8, ' ');
                xdmfFile << "<Attribute";
                xdmfFile << " Center=" << '"' << "Cell" << '"';
                xdmfFile << " Name=" << '"' << var_name << '"';
                xdmfFile << " AttributeType=" << '"' << "Scalar" << '"';
                xdmfFile << ">\n";
                xdmfFile << std::string(10, ' ');
                xdmfFile << "<DataItem";
                xdmfFile << " NumberType=" << '"' << "Float" << '"';
                xdmfFile << " Precision=" << '"' << precision << '"';

                xdmfFile << " Dimensions=\"";
                for (int idim = 2; idim >= 0; --idim) {
                    xdmfFile << ncells[idim];
                    xdmfFile << (idim == 0 ? "\"" : " ");
                }

                xdmfFile << " Format=" << '"' << "HDF" << '"';
                xdmfFile << ">\n";
                xdmfFile << std::string(12, ' ') << getFilename(prefix, it->first)
                         << ":/" << var_name << "\n";
                xdmfFile << std::string(10, ' ') << "</DataItem>\n";
                xdmfFile << std::string(8, ' ') << "</Attribute>\n";
            }
            // finalize grid file for the current time step
            xdmfFile << std::string(6, ' ') << "</Grid>\n";
        }

        // finalize Xdmf wrapper file
        xdmfFile << std::string(4, ' ') << "</Grid>\n";
        xdmfFile << std::string(2, ' ') << "</Domain>\n";
        xdmfFile << std::string(0, ' ') << "</Xdmf>\n";

        PDI_release("dl");
        PDI_release("origin");
        PDI_release("ncell");
        PDI_release("outputs_record");
        PDI_release("outputs_record_size");
    }
}

WriterPDI::WriterPDI(const UniformGrid& grid, const Params&,
                     const std::string& prefix,
                     const std::vector<std::pair<int, std::string>>&)
{
    std::array<int, 3> pdi_ncells;
    pdi_ncells[IX] = grid.m_nbCells[IX] * grid.m_dom[IX];
    pdi_ncells[IY] = grid.m_nbCells[IY] * grid.m_dom[IY];
    pdi_ncells[IZ] = grid.m_nbCells[IZ] * grid.m_dom[IZ];

    std::array<int, 3> pdi_ncells_local;
    pdi_ncells_local[IX] = grid.m_nbCells[IX];
    pdi_ncells_local[IY] = grid.m_nbCells[IY];
    pdi_ncells_local[IZ] = grid.m_nbCells[IZ];

    int tmp_rank=0;
#if defined(MPI_SESSION)
    MPI_Comm_rank(MPI_COMM_WORLD, &tmp_rank);
    m_mpi_coords = grid.comm.getCoords(grid.comm.rank());
#endif

    m_prefix=prefix;
    std::array<int, 3> pdi_start;
    pdi_start[IX] = grid.m_nbCells[IX] * m_mpi_coords[IX];
    pdi_start[IY] = grid.m_nbCells[IY] * m_mpi_coords[IY];
    pdi_start[IZ] = grid.m_nbCells[IZ] * m_mpi_coords[IZ];

    int prefix_size = prefix.size() + 1;
    int nvar = 9;

    std::array<Real, 3> origin;
    origin[IX] = grid.m_lowGlobal[IX];
    origin[IY] = grid.m_lowGlobal[IY];
    origin[IZ] = grid.m_lowGlobal[IZ];

    std::array<Real, 3> dl;
    dl[IX] = grid.m_dl[IX];
    dl[IY] = grid.m_dl[IY];
    dl[IZ] = grid.m_dl[IZ];

    int pdi_writer_time_step = 0;

    PDI_multi_expose("init_pdi_w_deisa",
                     "pdi_writer_time_step", &pdi_writer_time_step, PDI_OUT,
                     "mpi_coord", m_mpi_coords.data(), PDI_OUT,
                     "nvar", &nvar, PDI_OUT,
                     "ncell", pdi_ncells.data(), PDI_OUT,
                     "grid_size", pdi_ncells.data(), PDI_OUT,
                     "ghost", grid.m_ghostWidths.data(), PDI_OUT,
                     "ncell_local", pdi_ncells_local.data(), PDI_OUT,
                     "start", pdi_start.data(), PDI_OUT,
                     "origin", origin.data(), PDI_OUT,
                     "dl", dl.data(), PDI_OUT,
                     "restart_id", &m_restartId, PDI_OUT,
                     "prefix_size", &prefix_size, PDI_OUT,
                     "prefix", prefix.c_str(), PDI_OUT,
                     NULL);
    
}

void WriterPDI::write(HostConstArrayDyn u, const UniformGrid & grid,
                      Int iStep, Real time, Real gamma, Real mmw)
{
    std::array<int, 3> pdi_ncells;
    pdi_ncells[IX] = grid.m_nbCells[IX] * grid.m_dom[IX];
    pdi_ncells[IY] = grid.m_nbCells[IY] * grid.m_dom[IY];
    pdi_ncells[IZ] = grid.m_nbCells[IZ] * grid.m_dom[IZ];

    auto& outputId = WriterBase::m_outputId;

    char *prefix_c_str;
    PDI_access("prefix", (void **)&prefix_c_str, PDI_IN);
    std::string prefix(prefix_c_str);
    PDI_release("prefix");

    std::string filename = getFilename(prefix, outputId);
    int filename_size = filename.size();

    static int pdi_writer_time_step = 0;
    PDI_multi_expose("checkpoint",
                     "iStep", &iStep, PDI_OUT,
                     "pdi_writer_time_step", &pdi_writer_time_step, PDI_OUT,
                     "time", &time, PDI_OUT,
                     "Rstar_h", &code_units::constants::Rstar_h, PDI_OUT,
                     "gamma", &gamma, PDI_OUT,
                     "mmw", &mmw, PDI_OUT,
                     "output_id", &outputId, PDI_OUT,
                     "restart_id", &m_restartId, PDI_OUT,
                     "local_full_field", u.data(), PDI_OUT,
                     "filename_size", &filename_size, PDI_OUT,
                     "filename", filename.data(), PDI_OUT,
                     "grid_size", pdi_ncells.data(), PDI_OUT,
                     NULL);
    ++pdi_writer_time_step;

    WriterBase::m_previous_outputs.push_back(std::make_pair(outputId, time));

    ++outputId;

    int outputs_record_size = WriterBase::m_previous_outputs.size();

    PDI_multi_expose("write_xml",
                     "outputs_record_size", &outputs_record_size, PDI_OUT,
                     "outputs_record", WriterBase::m_previous_outputs.data(), PDI_OUT,
                     "restart_id", &m_restartId, PDI_OUT,
                     NULL);
}
}}
