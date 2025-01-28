#pragma once

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "ReaderBase.hpp"

#include <list>
#include <memory>
#include <utility>
#include <vector>

namespace hydro { namespace io
{

class Reader
{
public:
    Reader() = default;
    Reader(const UniformGrid& grid, const Params& params,
           const std::vector<std::pair<int, std::string>>& variables);
    Reader(const Reader& x) = default;
    Reader(Reader&& x) = default;
    ~Reader() = default;
    Reader& operator=(const Reader& x) = default;
    Reader& operator=(Reader&& x) = default;

    void read(HostArrayDyn u, const UniformGrid& grid,
              Int& iStep, Real& time, Int& outputId, Int& restartId) const;

private:
    std::shared_ptr<ReaderBase> io_reader = {};
};


}}
