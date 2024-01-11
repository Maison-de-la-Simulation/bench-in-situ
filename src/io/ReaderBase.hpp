#pragma once

#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

namespace hydro { namespace io
{

class ReaderBase
{
public:
    ReaderBase() = default;
    ReaderBase(const ReaderBase& x) = default;
    ReaderBase(ReaderBase&& x) = default;
    virtual ~ReaderBase() = default;
    ReaderBase& operator=(const ReaderBase& x) = default;
    ReaderBase& operator=(ReaderBase&& x) = default;

    virtual void read(HostArrayDyn u, const UniformGrid& grid,
                      Int& iStep, Real& time, Int& outputId, Int& restartId) = 0;
};

}}
