#include "WriterVTK.hpp"

#include "DistributedMemorySession.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "WriterBase.hpp"
#include "Utils.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace hydro { namespace io
{

WriterVTK::WriterVTK(const UniformGrid& grid, const Params& params,
                          const std::string& prefix,
                          const std::vector<std::pair<int, std::string>>& variables)
    : WriterBase {}
    , format      {params.output.format}
    , directory   {params.output.directory}
    , m_prefix    {prefix}
    , endian_type {utils::isLittleEndian() ? "LittleEndian" : "BigEndian"}
    , header_type {std::is_same<long unsigned, Int>::value ? "UInt64" : "UInt32"}
    , float_type  {std::is_same<double, Real>::value ? "Float64" : "Float32"}
    , xml_version {"1.0"}
    , vtk_version {"1.0"}
    , p_ext       {}
    , w_ext       {}
    , origin      {}
    , dl          {}
    , m_variables {variables}
{
#if defined(MPI_SESSION)
    std::array<int, three_d> coords {grid.comm.getCoords(grid.comm.rank())};
#else
    std::array<int, three_d> coords {};
#endif

    for (int idim=0; idim<3; ++idim)
    {
        p_ext[0][idim] = (idim<three_d ? grid.m_nbCells[idim]*coords[idim] : 0);
        p_ext[1][idim] = (idim<three_d ? grid.m_nbCells[idim]*(coords[idim]+1) : 1);
        w_ext[0][idim] = 0;
        w_ext[1][idim] = (idim<three_d ? grid.m_nbCells[idim]*params.mesh.dom[idim] : 1);
        origin[idim] = (idim<three_d ? grid.m_lowGlobal[idim] : constants::zero);
        dl[idim] = (idim<three_d ? grid.m_dl[idim] : constants::one);
    }
}

std::string WriterVTK::getFilename(Int outputId) const
{
    // write outputId in string outputNum
    std::ostringstream outputNum;
    outputNum << std::setw(std::numeric_limits<Int>::digits10);
    outputNum << std::setfill('0');
    outputNum << outputId;

#if defined(MPI_SESSION)
    std::string extension {"pvti"};
#else
    std::string extension {"vti"};
#endif

    // concatenate file prefix + file number + suffix
    std::string filename {m_prefix};
    filename += "_time_"+outputNum.str();
    filename += '.'+extension;
    return filename;
}

void WriterVTK::write(HostConstArrayDyn u, const UniformGrid& grid,
                           Int iStep, Real time, Real gamma, Real mmw)
{
    auto& outputId = WriterBase::m_outputId;

    WriterBase::m_previous_outputs.push_back(std::make_pair(outputId, time));
    if (Session::isIOProc())
    {
        writePvd();
#if defined(MPI_SESSION)
        writePvti(grid, outputId);
#endif
    }
    if (format == "ascii")
    {
        writeVtiAscii(u, grid, iStep, time, gamma, mmw);
    }
    else
    {
        writeVtiAppended(u, grid, iStep, time, gamma, mmw);
    }

    ++outputId;
}

void WriterVTK::writeVtiAscii(HostConstArrayDyn u, const UniformGrid& grid,
                                   Int iStep, Real time, Real gamma, Real mmw) const
{
    const auto& outputId = WriterBase::m_outputId;

    // open file
#if defined(MPI_SESSION)
    std::ofstream vtiFile {directory+'/'+getFilename(outputId, grid.comm.rank()), std::ofstream::trunc};
#else
    std::ofstream vtiFile {directory+'/'+getFilename(outputId), std::ofstream::trunc};
#endif
    if (!vtiFile.is_open())
    {
        std::cerr << "Error opening file, exiting..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    vtiFile << std::setprecision(std::numeric_limits<Real>::digits10);

    // write xml header
    vtiFile << "<?xml version=\"" << xml_version << "\"?>\n";

    // write VTKFile header
    vtiFile << "<VTKFile type=\"ImageData\"";
    vtiFile << " version="     << '"' << vtk_version << '"';
    vtiFile << " byte_order="  << '"' << endian_type << '"';
    vtiFile << " header_type=" << '"' << header_type << '"';
    vtiFile << ">\n";

    // write mesh extent
    vtiFile << std::string(2, ' ');
    vtiFile << "<ImageData";
    vtiFile << " WholeExtent=";
    vtiFile << '"' << p_ext[0][IX] << ' ' << p_ext[1][IX];
    vtiFile << ' ' << p_ext[0][IY] << ' ' << p_ext[1][IY];
    vtiFile << ' ' << p_ext[0][IZ] << ' ' << p_ext[1][IZ] << '"';
    vtiFile << " Origin=";
    vtiFile << '"' << origin[IX] << ' ' << origin[IY] << ' ' << origin[IZ] << '"';
    vtiFile << " Spacing=";
    vtiFile << '"' << dl[IX] << ' ' << dl[IY] << ' ' << dl[IZ] << '"';
    vtiFile << ">\n";

    vtiFile << std::string(4, ' ') << "<FieldData>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"gamma\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << gamma << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"mmw\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << mmw << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"Rstar_h\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << code_units::constants::Rstar_h << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"time\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << time << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"iStep\" NumberOfTuples=\"1\" type=\"Int64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << iStep << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"output_id\" NumberOfTuples=\"1\" type=\"Int64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << WriterBase::m_outputId << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"restart_id\" NumberOfTuples=\"1\" type=\"Int64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << WriterBase::m_restartId << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(4, ' ') << "</FieldData>\n";

    // write piece extent
    vtiFile << std::string(4, ' ');
    vtiFile << "<Piece";
    vtiFile << " Extent=";
    vtiFile << '"' << p_ext[0][IX] << ' ' << p_ext[1][IX];
    vtiFile << ' ' << p_ext[0][IY] << ' ' << p_ext[1][IY];
    vtiFile << ' ' << p_ext[0][IZ] << ' ' << p_ext[1][IZ] << '"';
    vtiFile << ">\n";

    // write point data
    vtiFile << std::string(6, ' ') << "<PointData>\n";
    vtiFile << std::string(6, ' ') << "</PointData>\n";

    // write cell data
    vtiFile << std::string(6, ' ');
    vtiFile << "<CellData";
    vtiFile << " Scalars=" << '"' << m_variables[0].second << '"';
    vtiFile << " Vectors=" << '"' << ""          << '"';
    vtiFile << ">\n";

    // write data array (ascii), remove ghost cells
    for (const auto& var : m_variables)
    {
        const int ivar = var.first;
        const std::string var_name = var.second;

        vtiFile << std::string(8, ' ') << "<DataArray";
        vtiFile << " Name="               << '"' << var_name   << '"';
        vtiFile << " NumberOfComponents=" << '"' << 1          << '"';
        vtiFile << " type="               << '"' << float_type << '"';
        vtiFile << " format="             << '"' << format     << '"';
        vtiFile << ">\n";

        vtiFile << std::string(9, ' ');
        for (Int index=0; index<grid.nbCells(); ++index)
        {
            if (grid.belongsToInnerDomain(index))
            {
                vtiFile << ' ';
                vtiFile << u(index, ivar);
            }
        }
        vtiFile << '\n';
        vtiFile << std::string(8, ' ') << "</DataArray>\n";
    } // end for ivar

    vtiFile << std::string(6, ' ') << "</CellData>\n";
    vtiFile << std::string(4, ' ') << "</Piece>\n";
    vtiFile << std::string(2, ' ') << "</ImageData>\n";
    vtiFile << std::string(0, ' ') << "</VTKFile>\n";
}

void WriterVTK::writeVtiAppended(HostConstArrayDyn u, const UniformGrid& grid,
                                      Int iStep, Real time, Real gamma, Real mmw) const
{
    const auto& outputId = WriterBase::m_outputId;

    // open file
#if defined(MPI_SESSION)
    std::ofstream vtiFile {directory+'/'+getFilename(outputId, grid.comm.rank()), std::ofstream::trunc};
#else
    std::ofstream vtiFile {directory+'/'+getFilename(outputId), std::ofstream::trunc};
#endif
    if (!vtiFile.is_open())
    {
        std::cerr << "Error opening file, exiting..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    vtiFile << std::setprecision(std::numeric_limits<Real>::digits10);

    // write VTKFile header
    vtiFile << "<VTKFile type=\"ImageData\"";
    vtiFile << " version="     << '"' << "1.0"       << '"';
    vtiFile << " byte_order="  << '"' << endian_type << '"';
    vtiFile << " header_type=" << '"' << header_type << '"';
    vtiFile << ">\n";

    // write mesh extent
    vtiFile << std::string(2, ' ');
    vtiFile << "<ImageData";
    vtiFile << " WholeExtent=";
    vtiFile << '"' << p_ext[0][IX] << ' ' << p_ext[1][IX];
    vtiFile << ' ' << p_ext[0][IY] << ' ' << p_ext[1][IY];
    vtiFile << ' ' << p_ext[0][IZ] << ' ' << p_ext[1][IZ] << '"';
    vtiFile << " Origin=";
    vtiFile << '"' << origin[IX] << ' ' << origin[IY] << ' ' << origin[IZ] << '"';
    vtiFile << " Spacing=";
    vtiFile << '"' << dl[IX] << ' ' << dl[IY] << ' ' << dl[IZ] << '"';
    vtiFile << ">\n";

    vtiFile << std::string(4, ' ') << "<FieldData>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"gamma\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << gamma << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"mmw\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << mmw << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"Rstar_h\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << code_units::constants::Rstar_h << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"time\" NumberOfTuples=\"1\" type=\"Float64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << time << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"iStep\" NumberOfTuples=\"1\" type=\"Int64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << iStep << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"output_id\" NumberOfTuples=\"1\" type=\"Int64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << WriterBase::m_outputId << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(6, ' ') << "<DataArray Name=\"restart_id\" NumberOfTuples=\"1\" type=\"Int64\" format=\"ascii\">\n";
    vtiFile << std::string(8, ' ') << WriterBase::m_restartId << '\n';
    vtiFile << std::string(6, ' ') << "</DataArray>\n";
    vtiFile << std::string(4, ' ') << "</FieldData>\n";

    // write piece extent
    vtiFile << std::string(4, ' ');
    vtiFile << "<Piece";
    vtiFile << " Extent=";
    vtiFile << '"' << p_ext[0][IX] << ' ' << p_ext[1][IX];
    vtiFile << ' ' << p_ext[0][IY] << ' ' << p_ext[1][IY];
    vtiFile << ' ' << p_ext[0][IZ] << ' ' << p_ext[1][IZ] << '"';
    vtiFile << ">\n";

    // write point data
    vtiFile << std::string(6, ' ') << "<PointData>\n";
    vtiFile << std::string(6, ' ') << "</PointData>\n";

    // write cell data
    Int inner_cells {1};
    for (int idim=0; idim<three_d; ++idim)
    {
        inner_cells *= grid.m_nbCells[idim];
    }

    vtiFile << std::string(6, ' ') << "<CellData>\n";
    for (const auto& var : m_variables)
    {
        const int ivar = var.first;
        const std::string var_name = var.second;

        const Int offset {ivar*(inner_cells*static_cast<Int>(sizeof(Real))+static_cast<Int>(sizeof(Int)))};
        vtiFile << std::string(8, ' ');
        vtiFile << "<DataArray";
        vtiFile << " Name="               << '"' << var_name   << '"';
        vtiFile << " NumberOfComponents=" << '"' << 1          << '"';
        vtiFile << " type="               << '"' << float_type << '"';
        vtiFile << " format="             << '"' << format     << '"';
        vtiFile << " offset="             << '"' << offset     << '"';
        vtiFile << "/>\n";
    } // end for ivar

    vtiFile << std::string(6, ' ') << "</CellData>\n";
    vtiFile << std::string(4, ' ') << "</Piece>\n";
    vtiFile << std::string(2, ' ') << "</ImageData>\n";

    vtiFile << std::string(2, ' ') << "<AppendedData encoding=\"raw\">\n";
    vtiFile << std::string(4, ' ') << '_';
    for (const auto& var : m_variables)
    {
        const int ivar = var.first;

        Int nbOfWords = inner_cells*static_cast<Int>(sizeof(Real));
        vtiFile.write(reinterpret_cast<char*>(&nbOfWords), sizeof(Int));
        for (Int j=0; j<grid.nbCells(); ++j)
        {
            if (grid.belongsToInnerDomain(j))
            {
                vtiFile.write(reinterpret_cast<char*>(&u(j, ivar)), sizeof(Real));
            }
        }
    }
    vtiFile << '\n';
    vtiFile << std::string(2, ' ') << "</AppendedData>\n";
    vtiFile << "</VTKFile>\n";
}

void WriterVTK::writePvd() const
{
    const auto& restartId = WriterBase::m_restartId;
    std::ostringstream restartNum;
    restartNum << std::setw(std::numeric_limits<Int>::digits10);
    restartNum << std::setfill('0');
    restartNum << restartId;

    std::string pvdFilename {directory+"/"+m_prefix+"_"+restartNum.str()+".pvd"};
    std::ofstream pvdFile    {pvdFilename, std::ofstream::trunc};
    if (!pvdFile.is_open())
    {
        std::cerr << "Error opening file, exiting..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    pvdFile << std::scientific;
    pvdFile << std::setprecision(std::numeric_limits<Real>::digits10);

    pvdFile << "<?xml version=\"" << xml_version << "\"?>\n";
    pvdFile << "<VTKFile type=" << '"' << "Collection" << '"';
    pvdFile << " version="      << '"' << vtk_version  << '"';
    pvdFile << " byte_order="   << '"' << endian_type  << '"';
    pvdFile << ">\n";

    pvdFile << std::string(2, ' ') << "<Collection>\n";
    for (const auto& previous_output : WriterBase::m_previous_outputs)
    {
        pvdFile << std::string(4, ' ');
        pvdFile << "<DataSet";
        pvdFile << " timestep=" << '"' << previous_output.second             << '"';
        pvdFile << " group="    << '"' << ""                     << '"';
        pvdFile << " part="     << '"' << 0                      << '"';
        pvdFile << " file="     << '"' << getFilename(previous_output.first) << '"';
        pvdFile << "/>\n";
    }
    pvdFile << std::string(2, ' ') << "</Collection>\n";

    pvdFile << "</VTKFile>\n";
}

#if defined(MPI_SESSION)
void WriterVTK::writePvti(const UniformGrid& grid, Int outputId) const
{
    // write outputId in string outputNum
    std::ostringstream outputNum;
    outputNum << std::setw(std::numeric_limits<Int>::digits10);
    outputNum << std::setfill('0');
    outputNum << outputId;

    // concatenate file prefix + file number + suffix
    std::string filename {directory+'/'+ m_prefix+"_time_"+outputNum.str()+".pvti"};

    // open file
    std::ofstream pvtiFile {filename, std::ofstream::trunc};
    if (!pvtiFile.is_open())
    {
        std::cerr << "Error opening file, exiting..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // write xml header
    pvtiFile << "<?xml version=\"" << xml_version << "\"?>\n";

    // write VTKFile header
    pvtiFile << "<VTKFile type=\"PImageData\"";
    pvtiFile << " version="     << '"' << vtk_version << '"';
    pvtiFile << " byte_order="  << '"' << endian_type << '"';
    pvtiFile << " header_type=" << '"' << header_type << '"';
    pvtiFile << ">\n";

    pvtiFile << std::string(2, ' ');
    pvtiFile << "<PImageData";
    pvtiFile << " WholeExtent=";
    pvtiFile << '"' << w_ext[0][IX] << ' ' << w_ext[1][IX];
    pvtiFile << ' ' << w_ext[0][IY] << ' ' << w_ext[1][IY];
    pvtiFile << ' ' << w_ext[0][IZ] << ' ' << w_ext[1][IZ] << '"';
    pvtiFile << " Origin=";
    pvtiFile << '"' << origin[IX] << ' ' << origin[IY] << ' ' << origin[IZ] << '"';
    pvtiFile << " Spacing=";
    pvtiFile << '"' << dl[IX] << ' ' << dl[IY] << ' ' << dl[IZ] << '"';
    pvtiFile << ">\n";

    pvtiFile << std::string(4, ' ');
    pvtiFile << "<PCellData Scalars=\"rho\">\n";
    for (const auto& var : m_variables)
    {
        const std::string var_name = var.second;

        pvtiFile << std::string(6, ' ');
        pvtiFile << "<PDataArray";
        pvtiFile << " type=\"" << float_type << '"';
        pvtiFile << " Name=\"" << var_name << '"';
        pvtiFile << " NumberOfComponents=\"" << 1 << "\"/>\n";
    }
    pvtiFile << std::string(4, ' ') << "</PCellData>\n";

    // one piece per MPI process
    for (int iPiece=0; iPiece<grid.comm.size(); ++iPiece)
    {
        // get MPI coords corresponding to MPI rank iPiece
        std::array<int, three_d> coords {grid.comm.getCoords(iPiece)};

        pvtiFile << std::string(4, ' ');
        pvtiFile << "<Piece Extent=";
        pvtiFile << '"' << coords[IX]*grid.m_nbCells[IX];
        pvtiFile << ' ' << (coords[IX]+1)*grid.m_nbCells[IX];
        pvtiFile << ' ' << (three_d>=2 ? coords[IY]*grid.m_nbCells[IY] : 0);
        pvtiFile << ' ' << (three_d>=2 ? (coords[IY]+1)*grid.m_nbCells[IY] : 1);
        pvtiFile << ' ' << (three_d==3 ? coords[IZ]*grid.m_nbCells[IZ] : 0);
        pvtiFile << ' ' << (three_d==3 ? (coords[IZ]+1)*grid.m_nbCells[IZ] : 1) << '"';
        pvtiFile << " Source=\"" << getFilename(outputId, iPiece) << "\"/>\n";
    }

    pvtiFile << std::string(2, ' ') << "</PImageData>\n";
    pvtiFile << std::string(0, ' ') << "</VTKFile>\n";
} // write_pvti_header

std::string WriterVTK::getFilename(Int outputId, int iPiece) const
{
    // write outputId in string outputNum
    std::ostringstream outputNum;
    outputNum << std::setw(std::numeric_limits<Int>::digits10);
    outputNum << std::setfill('0');
    outputNum << outputId;

    // write iPiece in string pieceNum
    std::ostringstream pieceNum;
    pieceNum << std::setw(std::numeric_limits<int>::digits10);
    pieceNum << std::setfill('0');
    pieceNum << iPiece;

    // concatenate file prefix + file number + suffix
    std::string filename {m_prefix};
    filename += "_time_"+outputNum.str();
    filename += "_piece_"+pieceNum.str();
    filename += ".vti";
    return filename;
}
#endif

class WriterVTK;

}}
