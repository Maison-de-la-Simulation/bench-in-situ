#pragma once

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "WriterBase.hpp"

#include <list>
#include <memory>
#include <vector>

namespace hydro { namespace io
{


struct WriterFactory
{
    static std::shared_ptr<WriterBase> New(const UniformGrid& grid, const Params& params,
                                                const std::string& type, const std::string& prefix,
                                                const std::vector<std::pair<int, std::string>>& variables);
};

struct WriterFactory;

}}
