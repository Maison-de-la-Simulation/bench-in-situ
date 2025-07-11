#include "WriterHDF5.hpp"
#include <hdf5.h> //??

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

WriterHDF5::WriterHDF5(const UniformGrid& grid, const Params&,
                        const std::string& prefix,
                        const std::vector<std::pair<int, std::string>>& variables)
    : WriterBase {}
    , m_prefix    {prefix}
    , hdf5_Int_type {std::is_same<long int, Int>::value ? H5T_NATIVE_INT64 : H5T_NATIVE_INT32}
    , hdf5_Real_type {std::is_same<double, Real>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT}
    , m_variables {variables}
{
    int tmp_rank=0;
#if defined(MPI_SESSION)
    MPI_Comm_rank(MPI_COMM_WORLD, &tmp_rank);
    m_mpi_coords = grid.comm.getCoords(grid.comm.rank());
#endif

    std::ostringstream mpi_prefix;
    mpi_prefix << std::setw(3) << std::setfill('0') << tmp_rank;    
    m_prefix.append("_r"+mpi_prefix.str());
}

std::string WriterHDF5::getFilename(Int outputId) const{

    // write outputId in string outputNum
    std::ostringstream outputNum;
    outputNum << std::setw(std::numeric_limits<Int>::digits10);
    outputNum << std::setfill('0');
    outputNum << outputId;

    // concatenate file prefix + file number + suffix
    std::string filename(m_prefix);
    filename += "_" + outputNum.str();
    filename += ".h5";
    return filename;
}

// write scalar_data
void WriterHDF5::write_scalar_data(const hid_t &file_id, const std::string & name, const hid_t &type, const void* data)
{
    raii_h5_hid const dataspace_id(H5Screate(H5S_SCALAR), H5Sclose);
    raii_h5_hid const dataset_id(H5Dcreate2(file_id, name.c_str(), type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);

    herr_t status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL,  H5P_DEFAULT, data);
    if (status < 0) {
        throw std::runtime_error("error in H5Dwrite for data "+name+"; status="+std::to_string(status));
    }
}

// write simple dataset(scalar, vector, array)
void WriterHDF5::write_simple_dataset_select(const hid_t &file_id, const std::string & name, const hid_t &type,
    const hsize_t dims_size, const hsize_t dims[/*size*/],
    const hsize_t hyperslab_start[/*size*/], const hsize_t hyperslab_count[/*size*/],
    const void* data)
{
    raii_h5_hid const dataspace_id(H5Screate_simple(dims_size, hyperslab_count, NULL), H5Sclose);
    raii_h5_hid const dataset_id(H5Dcreate2(file_id, name.c_str(), type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    raii_h5_hid const datamemory_id(H5Screate_simple(dims_size, dims, NULL), H5Sclose);

    /*
    * Define memory hyperslab.
    */
    herr_t status = H5Sselect_hyperslab(datamemory_id, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
    if (status < 0) {
        throw std::runtime_error("error in H5Sselect_hyperslab for data "+name+"; status="+std::to_string(status));
    }

    /*
    * Write the data to the dataset using default transfer properties.
    */
    status = H5Dwrite(dataset_id, type, datamemory_id, dataspace_id, H5P_DEFAULT, data);
    if (status < 0) {
        throw std::runtime_error("error in H5Dwrite for data "+name+"; status="+std::to_string(status));
    }
}

void WriterHDF5::writeXML(const UniformGrid & grid) const{
    if (!Session::isIOProc())
    {
        return;
    }
    const auto& restartId = WriterBase::m_restartId;
    std::ostringstream restartNum;
    restartNum << std::setw(std::numeric_limits<Int>::digits10);
    restartNum << std::setfill('0');
    restartNum << restartId;

    const std::string xdmfFilenameFull{m_prefix + '_' + restartNum.str() + ".xmf"};
    std::ofstream xdmfFile(xdmfFilenameFull, std::ofstream::trunc);
    if (!xdmfFile.is_open())
    {
        std::cerr << "Error opening file, exiting..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

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

    const std::vector< std::pair<Int, Real>> & outputs_record = WriterBase::m_previous_outputs;

    std::array<int, 3> ncells; // ncells of the global domain
    ncells[IX] = grid.m_nbCells[IX] * grid.m_dom[IX];
    ncells[IY] = grid.m_nbCells[IY] * grid.m_dom[IY];
    ncells[IZ] = grid.m_nbCells[IZ] * grid.m_dom[IZ];

    std::array<Real, 3> origin;
    origin[IX] = grid.m_lowGlobal[IX];
    origin[IY] = grid.m_lowGlobal[IY];
    origin[IZ] = grid.m_lowGlobal[IZ];

    std::array<Real, 3> dl;
    dl[IX] = grid.m_dl[IX];
    dl[IY] = grid.m_dl[IY];
    dl[IZ] = grid.m_dl[IZ];

    int precision = sizeof(Real);

    for (const auto & it: outputs_record) {
        xdmfFile << std::string(6, ' ');
        xdmfFile << "<Grid Name=" << '"' << "output" << '"';
        xdmfFile << " GridType=" << '"' << "Uniform" << '"';
        xdmfFile << ">\n";
        xdmfFile << std::string(8, ' ');
        xdmfFile << "<Time Value=" << '"' << it.second << '"'
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
        xdmfFile << std::string(12, ' ') << getFilename(it.first)
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
        xdmfFile << std::string(12, ' ') << getFilename(it.first)
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
        xdmfFile << std::string(12, ' ') << getFilename(it.first)
                    << ":/"
                    << "Rstar_h"
                    << "\n";
        xdmfFile << std::string(10, ' ') << "</DataItem>\n";
        xdmfFile << std::string(8, ' ') << "</Attribute>\n";

        for (const auto &var : m_variables) {
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
            xdmfFile << std::string(12, ' ') << getFilename(it.first)
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
}

void WriterHDF5::write(HostConstArrayDyn u, const UniformGrid& grid, 
                            Int iStep, Real time, Real gamma, Real mmw)
{
    auto& outputId = WriterBase::m_outputId;
    const auto& restartId = WriterBase::m_restartId;
    
    std::string filename = getFilename(outputId);
    const auto& Rstar_h = code_units::constants::Rstar_h;

    std::array<int, 3> ncells;
    ncells[IX] = grid.m_nbCells[IX] * grid.m_dom[IX];
    ncells[IY] = grid.m_nbCells[IY] * grid.m_dom[IY];
    ncells[IZ] = grid.m_nbCells[IZ] * grid.m_dom[IZ];

    const int dim = three_d;

    // create file 
    raii_h5_hid const file_id(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);

    //============================
    // write dimension of the grid_size
    std::vector<hsize_t> dims_grid(1);            // true dimension size of data in each direction
    std::vector<hsize_t> hyperslab_start_grid(1); // hyperslab_start in each direction
    dims_grid[0] = dim;
    hyperslab_start_grid[0] = 0;

    write_simple_dataset_select(file_id, "/grid_size", hdf5_Int_type, 1, dims_grid.data(), hyperslab_start_grid.data(), dims_grid.data(), ncells.data());

    //=========================
    // write physical variables
    //herr_t status=-1;

    // Create the data space for the dataset in memory and in file.
    std::vector<hsize_t> dims_var(dim);          // true dimension size of data in each direction
    std::vector<hsize_t> hyperslab_start(dim);   // hyperslab_start in each direction
    std::vector<hsize_t> dims_memory(dim);       // dimension size of the hyperslab in each direction (hyperslab_count)
    // The order of the dimension in hdf5 file is [grid.m_nbCells[2], grid.m_nbCells[1], grid.m_nbCells[0]]
    for (Int idim=0; idim<dim; ++idim) {
        dims_var[idim] = grid.m_nbCells[dim-1-idim] + 2*grid.m_ghostWidths[dim-1-idim];
        hyperslab_start[idim] = grid.m_ghostWidths[dim-1-idim];
        dims_memory[idim] = grid.m_nbCells[dim-1-idim];
    }

    Real *init_data;
    for (const auto& var : m_variables)
    {
        const int ivar = var.first;
        const std::string var_name = var.second;
        init_data=&u(0,ivar);
        // write var_name
        write_simple_dataset_select(file_id, var_name, hdf5_Real_type, dim, dims_var.data(), hyperslab_start.data(), dims_memory.data(), init_data);
    }
    init_data = nullptr;

    // write integer numbers
    write_scalar_data(file_id, "/iStep", hdf5_Int_type, &iStep);
    write_scalar_data(file_id, "/output_id", hdf5_Int_type, &outputId);
    write_scalar_data(file_id, "/restart_id", hdf5_Int_type, &restartId);

    // write real numbers
    write_scalar_data(file_id, "/Rstar_h", hdf5_Real_type, &Rstar_h);
    write_scalar_data(file_id, "/Time", hdf5_Real_type, &time);
    write_scalar_data(file_id, "/gamma", hdf5_Real_type, &gamma);
    write_scalar_data(file_id, "/mmw", hdf5_Real_type, &mmw);


    WriterBase::m_previous_outputs.push_back(std::make_pair(outputId, time));
    ++outputId;

    writeXML(grid);
}
    
}}
