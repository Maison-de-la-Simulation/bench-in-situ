#pragma once

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "ReaderBase.hpp"

#include <list>
#include <string>
#include <utility>
#include <vector>

namespace hydro { namespace io
{


class ReaderPDI : public ReaderBase
{
public:
    using RealVector    = RealVector3d;
    using IntVector     = IntVectorNd<three_d>;

    ReaderPDI() = default;
    ReaderPDI(const UniformGrid& grid, const Params& params,
                 const std::vector<std::pair<int, std::string>>& variables);
    ReaderPDI(const ReaderPDI& x) = default;
    ReaderPDI(ReaderPDI&& x) = default;
    ~ReaderPDI() override = default;

    void read(HostArrayDyn u, const UniformGrid& grid,
              Int& iStep, Real& time, Int& outputId, Int& restartId) final;

private:
    std::string filename;
    std::vector<std::pair<int, std::string>> m_variables;

}; // class ReaderPDI

}}
