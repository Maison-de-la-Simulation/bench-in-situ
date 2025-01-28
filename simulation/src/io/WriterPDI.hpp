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

    void write(HostConstArrayDyn u, const UniformGrid &grid,
               Int iStep, Real time, Real gamma, Real mmw) override;

    std::array<int, three_d> m_mpi_coords;
    std::string m_prefix;
};

}}
