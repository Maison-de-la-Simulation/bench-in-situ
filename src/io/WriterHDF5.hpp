#pragma once

#include <hdf5.h>

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "WriterBase.hpp"
#include "Utils.hpp"

#include <array>
#include <list>
#include <string>
#include <utility>
#include <vector>


namespace hydro { namespace io
{

class WriterHDF5 : public WriterBase
{
public:
    WriterHDF5() = default;
    WriterHDF5(const UniformGrid& grid, const Params& params,
              const std::string& prefix,
              const std::vector<std::pair<int, std::string>>& variables);
    WriterHDF5(const WriterHDF5& x) = default;
    WriterHDF5(WriterHDF5&& x) = default;
#if defined(__INTEL_COMPILER)
    ~WriterHDF5() {};
#else
    ~WriterHDF5() = default;
#endif // defined(__INTEL_COMPILER)
    //WriterHDF5& operator=(const WriterHDF5& x) = default;
    //WriterHDF5& operator=(WriterHDF5&& x) = default;

    void write(HostConstArrayDyn u, const UniformGrid& grid,
               Int iStep, Real time, Real gamma, Real mmw) override;

    std::string getFilename(Int outputId) const;

    private:
        void copy_data(std::vector<Real>& data, HostConstArrayDyn u, const UniformGrid & grid, Int ivar) const;
        void write_simple_dataset(const hid_t &file_id, const char* name, const hid_t &type, 
            const hsize_t size, const hsize_t dims[/*size*/], const void* data);
        void writeXML( const UniformGrid & grid) const;

        std::string m_prefix;
        hid_t hdf5_Int_type;
        hid_t hdf5_Real_type;
        
        std::vector<std::pair<int, std::string>> m_variables; // list of physical array write in the file
        std::array<int, three_d> m_mpi_coords;
};

class WriterHDF5;
}}
