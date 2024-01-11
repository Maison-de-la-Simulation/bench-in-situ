#pragma once

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

class WriterPDI : public WriterBase
{
public:
    WriterPDI() = default;
    WriterPDI(const UniformGrid& grid, const Params& params,
              const std::string& prefix,
              const std::vector<std::pair<int, std::string>>& variables);
    WriterPDI(const WriterPDI& x) = default;
    WriterPDI(WriterPDI&& x) = default;
    ~WriterPDI() override = default;
    //WriterPDI& operator=(const WriterPDI& x) = default;
    //WriterPDI& operator=(WriterPDI&& x) = default;

    void write(HostConstArrayDyn u, const UniformGrid& grid,
               Int iStep, Real time, Real gamma, Real mmw) override;
    void write_mean(global_mean means) override;
    void write_profile(void* profile_data, const UniformGrid& grid) override;
    void write_slice(const UniformGrid& grid, void* h_slice_data, void* v_slice_data, bool contains_middle_z) final;

    std::array<int, three_d> m_mpi_coords;
    std::string m_prefix;
};

}}
