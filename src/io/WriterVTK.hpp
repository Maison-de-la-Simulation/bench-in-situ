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

class WriterVTK : public WriterBase
{
public:
    WriterVTK() = default;
    WriterVTK(const UniformGrid& grid, const Params& params,
              const std::string& prefix,
              const std::vector<std::pair<int, std::string>>& variables);
    WriterVTK(const WriterVTK& x) = default;
    WriterVTK(WriterVTK&& x) = default;
#if defined(__INTEL_COMPILER)
    ~WriterVTK() {};
#else
    ~WriterVTK() = default;
#endif // defined(__INTEL_COMPILER)
    //WriterVTK& operator=(const WriterVTK& x) = default;
    //WriterVTK& operator=(WriterVTK&& x) = default;

    void write(HostConstArrayDyn u, const UniformGrid& grid,
               Int iStep, Real time, Real gamma, Real mmw) final;
    void writeVtiAscii(HostConstArrayDyn u, const UniformGrid& grid,
                       Int iStep, Real time, Real gamma, Real mmw) const;
    void writeVtiAppended(HostConstArrayDyn u, const UniformGrid& grid,
                          Int iStep, Real time, Real gamma, Real mmw) const;
    void writePvd() const;
    std::string getFilename(Int iStep) const;
#if defined(MPI_SESSION)
    void writePvti(const UniformGrid& grid, Int iStep) const;
    std::string getFilename(Int iStep, int iPiece) const;
#endif

private:
    std::string format;
    std::string directory;
    std::string m_prefix;
    std::string endian_type;
    std::string header_type;
    std::string float_type;
    std::string xml_version;
    std::string vtk_version;
    std::array<std::array<Int, 3>, 2> p_ext;
    std::array<std::array<Int, 3>, 2> w_ext;
    std::array<Real, 3> origin;
    std::array<Real, 3> dl;
    std::vector<std::pair<int, std::string>> m_variables;
};
class WriterVTK;
}}
