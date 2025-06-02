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

// write simple dataset(scalar, vector, array)
void WriterHDF5::write_simple_dataset_select(const hid_t &file_id, const char* name, const hid_t &type,
    const hsize_t dims_size, const hsize_t dims[/*size*/],
    const hsize_t hyperslab_start[/*size*/], const hsize_t hyperslab_count[/*size*/],
    const void* data)
{
    // If we write all element in data,  hyperslab_start = 0  and hyperslab_count = dims
    herr_t status=-1;
    herr_t status_tmp=0;
    hid_t dataset_id;
    hid_t dataspace_id;
    hid_t datamemory_id;

    if (dims_size == 1 && dims[0] == 1){
        dataspace_id = H5Screate(H5S_SCALAR);
    }
    else{
        /*
        * Describe the size of the array and create the data space for fixed
        * size dataset.
        */
        dataspace_id = H5Screate_simple(dims_size, hyperslab_count, NULL);
    }
    if (dataspace_id < 0) {
        std::cout << "error "<< dataspace_id <<" in creating dataspace for " << name << std::endl;
        goto end_function;
    }

    /*
    * Create a new dataset within the file using defined dataspace and
    * datatype and default dataset creation properties.
    */
    dataset_id = H5Dcreate2(file_id, name, type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cout << "error in H5Dcreate2 for "<< name <<"; dataset_id=" << dataset_id << std::endl;
        goto close_dataspace_id;
    }

    datamemory_id = H5Screate_simple(dims_size, dims, NULL);
    if (datamemory_id < 0) {
        std::cout << "error in H5Dcreate2 for "<< name <<"; datamemory_id=" << datamemory_id << std::endl;
        goto close_dataset_id;
    }

    /*
    * Define memory hyperslab.
    */
    status = H5Sselect_hyperslab(datamemory_id, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
    if (status < 0) {
        std::cout << "error in H5Sselect_hyperslab for "<< name <<"; status=" << status << std::endl;
        goto close_datamemory_id;
    }

    /*
    * Write the data to the dataset using default transfer properties.
    */
    status = H5Dwrite(dataset_id, type, datamemory_id, dataspace_id, H5P_DEFAULT, data);
    if (status < 0) {
        std::cout << "error in H5Dwrite for "<< name <<"; status=" << status << std::endl;
    }

close_datamemory_id:
    status_tmp=H5Sclose(datamemory_id);
    if (status_tmp < 0) {
        status=status_tmp;
        std::cout << "error in closing dataset for "<< name <<"; status=" << status_tmp << std::endl;
    }

close_dataset_id:
    status_tmp=H5Dclose(dataset_id);
    if (status_tmp < 0) {
        status=status_tmp;
        std::cout << "error in closing dataset for "<< name <<"; status=" << status_tmp << std::endl;
    }

close_dataspace_id:
    status_tmp=H5Sclose(dataspace_id);
    if (status_tmp < 0) {
        status=status_tmp;
        std::cout << "error in closing dataspace for "<< name <<"; status=" << status_tmp << std::endl;
    }

end_function:
    if (status < 0) {
        exit(EXIT_FAILURE);
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

    // Jacques:::: A voir si on met la version local
    std::array<int, 3> ncells;
    ncells[IX] = grid.m_nbCells[IX] * grid.m_dom[IX];
    ncells[IY] = grid.m_nbCells[IY] * grid.m_dom[IY];
    ncells[IZ] = grid.m_nbCells[IZ] * grid.m_dom[IZ];

    const int dim = three_d;

    // create file 
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id < 0) {
        std::cout << "Error "<< file_id << " in H5Fcreate for filename=" << filename << std::endl;
		exit(1);
	}

    //============================
    // write dimension of the grid_size
    
    hsize_t dims_grid[1];
    dims_grid[0] = dim;
    hsize_t hyperslab_start_grid[1];
    hyperslab_start_grid[0] = 0;

    write_simple_dataset_select(file_id, "/grid_size", hdf5_Int_type, 1, dims_grid, hyperslab_start_grid, dims_grid, ncells.data());

    //=========================
    // write physical variables
    herr_t status=-1;

    // Create the data space for the dataset in memory and in file.
    hsize_t dims_var[dim];          // true dimension size of data in each direction
    hsize_t hyperslab_start[dim];   // hyperslab_start in each direction
    hsize_t dims_memory[dim];       // dimension size of the hyperslab in each direction (hyperslab_count)
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
        write_simple_dataset_select(file_id, var_name.c_str(), hdf5_Real_type, dim, dims_var, hyperslab_start, dims_memory, init_data);
    }
    init_data = nullptr;

    hsize_t dims_scalar[1];
    dims_scalar[0] = 1;
    hsize_t hyperslab_start_scalar[1];
    hyperslab_start_scalar[0] = 0;
    
    // write integer numbers
    write_simple_dataset_select(file_id, "/iStep", hdf5_Int_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &iStep);
    write_simple_dataset_select(file_id, "/output_id", hdf5_Int_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &outputId);
    write_simple_dataset_select(file_id, "/restart_id", hdf5_Int_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &restartId);
   
    // write real numbers
    write_simple_dataset_select(file_id, "/Rstar_h", hdf5_Real_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &Rstar_h);
    write_simple_dataset_select(file_id, "/Time", hdf5_Real_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &time);
    write_simple_dataset_select(file_id, "/gamma", hdf5_Real_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &gamma);
    write_simple_dataset_select(file_id, "/mmw", hdf5_Real_type, 1, dims_scalar, hyperslab_start_scalar, dims_scalar, &mmw);

    // close file
    status = H5Fclose(file_id);
	if (status < 0) {
        std::cout << "Error "<< status << " in H5Fclose for filename=" << filename << std::endl;
		exit(1);
	}

    WriterBase::m_previous_outputs.push_back(std::make_pair(outputId, time));
    ++outputId;

    writeXML(grid);
}
    
}}
