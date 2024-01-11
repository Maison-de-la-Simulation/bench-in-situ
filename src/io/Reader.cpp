#include "Reader.hpp"

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "ReaderBase.hpp"
#if defined(Euler_ENABLE_PDI)
#include "ReaderPDI.hpp"
#endif

#include <list>
#include <memory>
#include <utility>
#include <vector>

namespace hydro { namespace io
{

#if defined(Euler_ENABLE_PDI)

Reader::Reader(const UniformGrid& grid, const Params& params,
                    const std::vector<std::pair<int, std::string>>& variables)
{
    io_reader = std::make_shared<ReaderPDI>(grid, params, variables);
}


void Reader::read(HostArrayDyn u, const UniformGrid& grid,
                       Int& iStep, Real& time, Int& outputId, Int& restartId) const
{
    io_reader->read(u, grid, iStep, time, outputId, restartId);
}

#else

Reader::Reader(const UniformGrid&, const Params&,
                    const std::vector<std::pair<int, std::string>>&)
{
}


void Reader::read(HostArrayDyn, const UniformGrid&,
                       Int&, Real&, Int&, Int&) const
{
}

#endif

}}
