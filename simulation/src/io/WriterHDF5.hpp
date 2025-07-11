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
#include <functional>


namespace hydro { namespace io
{

/** A RAII-style wrapper for HDF5 hid_t.
 *
 * This calls the provided destroyer function when the hid_t goes out of scope.
 */
class raii_h5_hid
{
private:
    /// The wrapped hid_t
    hid_t m_id;

    /// The destroyer function rto call, or null if none
    std::function<herr_t(hid_t)> m_destroyer;

public:
    raii_h5_hid(hid_t id, herr_t (*f)(hid_t)) : m_id(id), m_destroyer(f)
    {
        if (m_id < 0 || !m_destroyer) {
            throw std::runtime_error("bench-in-situ error: creating h5 id failed");
        }
    }

    /// No copy possible
    raii_h5_hid(const raii_h5_hid&) = delete;
    raii_h5_hid(raii_h5_hid&&) = delete;

    ~raii_h5_hid() noexcept
    {
        if (m_id >= 0 && m_destroyer) {
            try {
                m_destroyer(m_id);
            }
            catch(...){
                std::cerr << "bench-in-situ error: closing raii_h5_hid failed" << std::endl;
            }
        }
    }

    /// No copy possible
    raii_h5_hid& operator=(const raii_h5_hid&) = delete;
    raii_h5_hid& operator=(raii_h5_hid&&) = delete;


    /** Supports using the Raii_5d_hid as a raw hid_t.
	 * \return the raw hid_t
	 */
	operator hid_t () const noexcept
    {
        return m_id;
    }
};

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
        void write_scalar_data(const hid_t &file_id, const std::string & name, const hid_t &type, const void* data);
        void write_simple_dataset_select(const hid_t &file_id, const std::string& name, const hid_t &type,
            const hsize_t dims_size, const hsize_t dims[/*size*/], const hsize_t hyperslab_start[/*size*/],
            const hsize_t hyperslab_count[/*size*/], const void* data);
        void writeXML( const UniformGrid & grid) const;

        std::string m_prefix;
        hid_t hdf5_Int_type;
        hid_t hdf5_Real_type;
        
        std::vector<std::pair<int, std::string>> m_variables; // list of physical array write in the file
        std::array<int, three_d> m_mpi_coords;
};

class WriterHDF5;
}}
