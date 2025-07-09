#pragma once

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "WriterBase.hpp"
#include "Utils.hpp"
#include "Timer.hpp"

#include <array>
#include <list>
#include <string>
#include <utility>
#include <vector>



namespace hydro { namespace io
{

// std::string getFilename(std::string const &prefix, Int outputId);

class WriterGpuPDI : public WriterBase
{
public:
    WriterGpuPDI() = default;
    WriterGpuPDI(const UniformGrid& grid, const Params& params,
              const std::string& prefix,
              const std::vector<std::pair<int, std::string>>& variables);
    WriterGpuPDI(const WriterGpuPDI& x) = default;
    WriterGpuPDI(WriterGpuPDI&& x) = default;
    ~WriterGpuPDI() override = default;
    //WriterGpuPDI& operator=(const WriterGpuPDI& x) = default;
    //WriterGpuPDI& operator=(WriterGpuPDI&& x) = default;

    void write(HostConstArrayDyn u, const UniformGrid &grid,
               Int iStep, Real time, Real gamma, Real mmw) override;

    std::array<int, three_d> m_mpi_coords;
    std::string m_prefix;

private:
    DebugTimer debugTimer;
};

}}
